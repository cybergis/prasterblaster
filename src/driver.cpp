/*!
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE 
 *
 * This software is in the public domain, furnished "as is", without
 * technical support, and with no warranty, express or implied, as to
 * its usefulness for any purpose.
 *
 * @section DESCRIPTION
 *
 * 
 */

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cctype>

#include <mpi.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mollweide.h"
#include "projectedraster.hh"
#include "rasterchunk.hh"
#include "reprojector.hh"

#include <gdal_priv.h>
#include <ogr_spatialref.h>

using std::shared_ptr;

int driver(string input_raster, 
	   string output_filename, 
	   string output_srs, 
	   string fillvalue,
	   int partition_count)
{
	int rank = 0;
	int process_count = 1;
	shared_ptr<ProjectedRaster> in, out;
	shared_ptr<Projection> in_proj, out_proj;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);

	// Open input raster and check for errors
	if (rank == 0) {
	  printf("Opening input raster...");
	  fflush(stdout);
	}
        in = shared_ptr<ProjectedRaster>(new ProjectedRaster(input_raster));
	if (rank == 0) 
	  printf("done\n");

        if (in->isReady() == true) {

        } else {
                fprintf(stderr, "Error opening input raster\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
                return 1;
        }

	// If we are rank 0, create the output projection and raster
	if (rank == 0) {
		in_proj = shared_ptr<Projection>(in->getProjection());
		OGRSpatialReference srs; 

		OGRErr err = srs.importFromProj4(output_srs.c_str());

		if (err != OGRERR_NONE) {
			fprintf(stderr, "Error parsing projection!\n");
			return -1;
		}
		long proj_code, datum_code, zone;
		double *params = NULL;
		
		srs.exportToUSGS(&proj_code, &zone, &params, &datum_code);
		
		out_proj = shared_ptr<Projection>(Transformer::convertProjection(SNSOID));

		out_proj->setUnits(in_proj->units());
		out_proj->setDatum(in_proj->datum());
		out_proj->setParams(params);

		OGRFree(params);

		printf("Creating output raster...");
		fflush(stdout);
		bool result  = ProjectedRaster::CreateRaster(output_filename,
							     in,
							     shared_ptr<Projection>(out_proj->copy()),
							     in->type ,
							     in->pixel_size);

		
		if (result == false) {
			fprintf(stderr, "Failed to create output raster!\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		} else {
			printf("done\n");
		}
		printf("Syncing nodes...");
		MPI_Barrier(MPI_COMM_WORLD);
		printf("done\n");
	} else {
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	
	// Now we re-open the output raster on each node.
	if (rank == 0) {
		printf("Opening new output raster...");
		fflush(stdout);
	}
	out = shared_ptr<ProjectedRaster>(new ProjectedRaster(output_filename));
	if (out == 0) {
		fprintf(stderr, "Output allocation failed, something is very wrong!\n");
		return 1;

	}
	
	if (!in->isReady()) {
		fprintf(stderr, "Error opening input raster, not ready!\n");
		return 1;

	}

	if (!out->isReady()) {
		fprintf(stderr, "Error opening output raster, not ready!\n");
		return 1;

	}
	if (rank == 0) 
		printf("done\n");

	// Now preform actual reprojection
	if (rank == 0) {
		printf("Reprojecting...");
		fflush(stdout);
	}

	std::vector<Area> part_areas = PartitionByCount(out, partition_count);
	RasterChunk::RasterChunk *out_chunk, *in_chunk;
	Area output_area, input_area;
	int first_index, last_index;

	first_index = rank * (part_areas.size()/process_count);
	last_index = ((rank+1) * (part_areas.size()/process_count)) - 1;

	// If we are the last process, make sure all of the
	// RasterChunks are allocated
	if (rank == process_count -1) {
		last_index = part_areas.size() - 1;
	}


	extern int analyze_partitions;
	if (analyze_partitions == 1) {
		for (int i = 0; i < partition_count; ++i) {
			input_area = MapDestinationAreatoSource(out, in, part_areas.at(i));
			Area part = part_areas.at(i);
			out_chunk = out->createEmptyRasterChunk(part_areas.at(i));
			in_chunk = in->createEmptyRasterChunk(input_area);
			printf("\n-------------------------------------------------------------------------------\n");
			printf("OUTPUT AREA: UL: %f %f LR: %f %f \n", part.ul.x, part.ul.x, part.lr.x, part.lr.y);
			printf("             --- %d rows %d columns\n", out_chunk->row_count_, out_chunk->column_count_);
			printf("INPUT AREA: UL: %f %f LR: %f %f \n", input_area.ul.x, input_area.ul.y, input_area.lr.x, input_area.lr.y);
			printf("             --- %d rows %d columns\n", in_chunk->row_count_, in_chunk->column_count_);
			printf("\n-------------------------------------------------------------------------------\n");
			delete out_chunk;
			delete in_chunk;
			       
		}
		MPI_Abort(MPI_COMM_WORLD, 0);
	}

	for (int i = first_index; i <= last_index; ++i) {
		output_area = part_areas.at(i);
		input_area = MapDestinationAreatoSource(out, in, output_area);
		out_chunk = out->createAllocatedRasterChunk(output_area);
		in_chunk = in->createRasterChunk(input_area);

		if (in_chunk == NULL) { // break 
			fprintf(stderr, "Input RasterChunk allocation error!\n");
			return 1;
		}

		if (out_chunk == NULL) { // break 
			fprintf(stderr, "Output RasterChunk allocation error!\n");
			return 1;
		}
		printf("%d,", i);
		fflush(stdout);
		if (ReprojectChunk(in_chunk, out_chunk) == false) {
			fprintf(stderr, "Error reprojecting chunk #%d\n", i);
		}

		// Now write RasterChunk to output
		if (out->writeRasterChunk(out_chunk) == false) {
			fprintf(stderr, "Error writing chunk!\n");
		} 

		

		// Cleanup
		delete out_chunk;
		delete in_chunk;
		out_chunk = NULL;
		in_chunk = NULL;

	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
		printf(" done!\n");

	// Cleanup

	return 0;
}

