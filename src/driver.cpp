
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
#include <sys/types.h>
#include <sys/stat.h>

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
	char *output_template = strdup(temporary_path.c_str());
	std::vector<string> chunk_names;
	PRB_ERROR result;

	long out_projcode, out_datumcode;
	long out_datum;
	double *out_params;

	std::vector<std::string> temp_output_files;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);

	double begin_prologue = MPI_Wtime();

	std::stringstream sout;
	sout << rank;
	// Create directory for this processor's temp files
	temporary_path = temporary_path + string("/") + sout.str();
	mkdir(temporary_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
	
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

	// Create an output raster for rank 0
	if (rank == 0) {
		printf("Creating output raster...");
		fflush(stdout);
		result = CreateOutputRaster(in, output_filename, in->pixel_size, output_srs);
	}

	if (result != NO_ERROR) {
		fprintf(stderr, "Failed to create output raster!\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	} else {
		out.reset(new ProjectedRaster(output_filename));
		out_proj = out->getProjection();
		printf("done\n");
	}

	if (!in->isReady()) {
		fprintf(stderr, "Error opening input raster, not ready!\n");
		return 1;
	}

	// Now perform actual reprojection
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

	double end_prologue = MPI_Wtime();

	for (int i = first_index; i <= last_index; ++i) {
		
		// Create temporary file for output chunk
		string out_path = temporary_path + "/prbXXXXXXXX";
		char *tempfilename = strdup(out_path.c_str());
		int fd = mkstemp(tempfilename);
		out_path = tempfilename;
		temp_output_files.push_back(out_path);

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
			// If we get this error it's probably a bug in RasterMinbox
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
		if (WriteRasterChunk(out_path, out_chunk) != NO_ERROR) {
			fprintf(stderr, "Rank %d, Error writing chunk!\n", rank);
		} 

		

		// Cleanup
		delete out_chunk;
		delete in_chunk;
		out_chunk = NULL;
		in_chunk = NULL;

	}

	double end_reprojection = MPI_Wtime();


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


	double end_epilogue = MPI_Wtime();

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		double pro_time = end_prologue - begin_prologue;
		double reproject_time = end_reprojection - end_prologue;
		double epilogue_time = end_epilogue - end_reprojection;
		double total_time = pro_time + reproject_time + epilogue_time;
		printf("\n\n\n");
		printf("Reprojection Complete! ***Run times** |Prologue: %f| Reprojection %f| Epilogue: %f| TOTAL %f SECONDS",
		       pro_time, reproject_time, epilogue_time, total_time);
		printf("\n\n\n");
	}
	return 0;
}

