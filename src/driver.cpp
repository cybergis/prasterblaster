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
#include <sstream>
#include <unistd.h>

#include <mpi.h>

#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mollweide.h"
#include "projectedraster.hh"
#include "rasterchunk.hh"
#include "reprojector.hh"
#include "sharedptr.hh"


			 

int driver(string input_raster, 
	   string output_filename, 
	   string temporary_path,
	   string output_srs, 
	   string resampler,
	   string fillvalue,
	   int partition_count)
{
	int rank = 0;
	int process_count = 1;
	shared_ptr<ProjectedRaster> in, out;
	shared_ptr<Projection> in_proj, out_proj;
	string final_output_filename = output_filename;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);

	std::stringstream sout;
	sout << rank;
	if (rank != 0) {
		output_filename = output_filename + sout.str();
	}
	
	// Open Input raster and check for errors
	if (rank == 0) {
	  printf("Opening input raster...");
	  fflush(stdout);
	}
        in = shared_ptr<ProjectedRaster>(new ProjectedRaster(input_raster));
	if (rank == 0) 
	  printf("done\n");

        if (in->isReady() == false) {
                fprintf(stderr, "Error opening input raster\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
                return 1;
        }

	// Create an output raster for each process
	bool result = CreateOutputRaster(in, output_filename, output_srs);
	if (result == false) {
		fprintf(stderr, "Failed to create output raster!\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	} else {
		printf("done\n");
	}
	printf("Syncing nodes...");
	MPI_Barrier(MPI_COMM_WORLD);
	printf("done\n");
	
	
	
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
	
	if (rank == 0) {
		printf("Finding partitions...");
		fflush(stdout);
	}
	std::vector<Area> part_areas = PartitionByCount(out, partition_count);
	if (rank == 0) {
		printf("done!\n");
	}
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


	int analyze_partitions = 0;
	if (analyze_partitions == 1) {
		for (int i = 0; i < partition_count; ++i) {
			input_area = RasterMinbox(out, in, part_areas.at(i));
			Area part = part_areas.at(i);
			out_chunk = out->createEmptyRasterChunk(part_areas.at(i));
			in_chunk = in->createEmptyRasterChunk(input_area);
			printf("\n-------------------------------------------------------------------------------\n");
			printf("PARTITION #%d\n", i);
			printf("OUTPUT AREA: UL: %f %f LR: %f %f \n", part.ul.x, part.ul.y, part.lr.x, part.lr.y);
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

		// Swap y dimension
		output_area.ul.y = out->getRowCount() - output_area.ul.y - 1;
		output_area.lr.y = out->getRowCount() - output_area.lr.y - 1;
		part_areas.at(i) = output_area;
		if (output_area.ul.x == -1) {
			continue;
		}
		
		input_area = RasterMinbox(out, in, output_area);
		out_chunk = out->createAllocatedRasterChunk(output_area);
		in_chunk = in->createRasterChunk(input_area);

		if (in_chunk == NULL) { // break 
			fprintf(stderr, "Rank: %d: Input RasterChunk allocation error!\n", rank);
			fprintf(stderr, "      %f %f %f %f\n", input_area.ul.x, input_area.ul.y,
				input_area.lr.x, input_area.lr.y);
			return 1;
		}

		if (out_chunk == NULL) { // break 
			fprintf(stderr, "Rank %d: Output RasterChunk allocation error!\n", rank);
			return 1;
		}
		int t = last_index - first_index;
		if ((t/100) == 0 || i % (t/100) == 0) {
			int c = i - first_index;
			for (int k=0; k<rank; ++k) {
				printf("\t\t\t\t");
			}
			printf("[Rank %d] %d/%d {%.2f%%}\n", 
			       rank, c, t, 100*((double)c/(double)t));
			fflush(stdout);
		}
		if (ReprojectChunk(in_chunk, out_chunk, fillvalue, resampler) == false) {
			fprintf(stderr, "Rank %d, Error reprojecting chunk #%d\n", rank, i);
		}

		// Now write RasterChunk to output
		if (out->writeRasterChunk(out_chunk) == false) {
			fprintf(stderr, "Rank %d, Error writing chunk!\n", rank);
		} 

		

		// Cleanup
		delete out_chunk;
		delete in_chunk;
		out_chunk = NULL;
		in_chunk = NULL;

	}

	// Now copy temporary rasters to rank 0's 
	for (int i=1; i<process_count; ++i) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == i) {
			printf("Copying from rank %d...", rank);
			fflush(stdout);
			shared_ptr<ProjectedRaster> main_raster(new ProjectedRaster(final_output_filename));
			for (int i = first_index; i <= last_index; ++i) {
				out_chunk = out->createRasterChunk(part_areas.at(i));
				if (out_chunk == NULL) {
					fprintf(stderr, "Error Allocating output chunk!\n");
					exit(1);
				}
				main_raster->writeRasterChunk(out_chunk);
				delete out_chunk;
				unlink(output_filename.c_str());
			}			
			printf("done!\n");
			main_raster.reset();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
		printf(" done!\n");

	// Cleanup

	return 0;
}

