/*
  Programmer: David Mattli
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

using std::shared_ptr;

double params[15] =  { 6370997.000000, 
		       0, 0, 0, 0, 
		       0, 0, 0, 0, 0, 
		       0, 0, 0, 0, 0};


int driver(string input_raster, string output_filename, string output_projection)
{
	int rank = 0;
	int process_count = 1;
	shared_ptr<ProjectedRaster> in, out;
	shared_ptr<Projection> in_proj, out_proj;
	int partition_count = 1;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);
	partition_count = process_count * 2;

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
		// Downcase output projection name
		std::transform(output_projection.begin(),
			       output_projection.end(),
			       output_projection.begin(),
			       (int(*) (int))std::tolower);

		in_proj = shared_ptr<Projection>(in->getProjection());

		if (output_projection == "hammer") {
			out_proj = shared_ptr<Projection>(Transformer::convertProjection(HAMMER));
		} else if (output_projection == "mollweide") {
			out_proj = shared_ptr<Projection>(Transformer::convertProjection(MOLL));
		} else if (output_projection == "sinosoidal") {
			return 1;
		} else {
			// fail...
			return 1;
		}
		out_proj->setUnits(in_proj->units());
		out_proj->setDatum(in_proj->datum());
		out_proj->setParams(in_proj->params());

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
		
		printf("Reprojecting UL: %f %f rows %d cols %d\n", in_chunk->raster_location_.x, in_chunk->raster_location_.y,
		       in_chunk->row_count_, in_chunk->column_count_);

		if (ReprojectChunk(*in_chunk, *out_chunk) == false) {
			fprintf(stderr, "Error reprojecting chunk #%d\n", i);
		}

		// Now write RasterChunk to output
		if (out->writeRasterChunk(out_chunk) == false) {
			fprintf(stderr, "Error writing chunk!\n");
		} 

		

		// Cleanup
		delete out_chunk;
		delete in_chunk;
	}

	if (rank == 0)
		printf("done!\n");

	// Cleanup

	return 0;
}

