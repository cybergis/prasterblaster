//
// Copyright 0000 <Nobody>
// @file
// @author David Matthew Mattli <dmattli@usgs.gov>
//
// @section LICENSE
//
// This software is in the public domain, furnished "as is", without
// technical support, and with no warranty, express or implied, as to
// its usefulness for any purpose.
//
// @section DESCRIPTION
//
// This file demonstrates how to use librasterblaster to implement parallel
// raster reprojection. This implementation uses a MPI I/O to write to a tiff
// file in parallel.
//
//

#include <vector>

#include "configuration.h"
#include "projectedraster.h"
#include "quadtree.h"
#include "reprojection_tools.h"
#include "sharedptr.h"

#include "sptw.h"

using namespace librasterblaster;
using namespace sptw;

int main(int argc, char *argv[]) {
  int rank = 0;
  int process_count = 0;
  RasterChunk *in_chunk, *out_chunk;
  
  // Give MPI_Init first run at the command-line arguments
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);
  
  // Initialize Configuration object
  Configuration conf(argc, argv);
  
  if (conf.input_filename == "" || conf.output_filename == "") {
    fprintf(stderr, "Specify an input and output filename\n");
    return 0;
  }

  if (conf.partition_size == 0) {
    conf.partition_size = 50000;
  }
  
  // Open the input raster
  shared_ptr<ProjectedRaster> input_raster(new ProjectedRaster(conf.input_filename));
  if (input_raster->ready() == false) {
    fprintf(stderr, "Error opening input raster!\n");
    return 0;
  }

  if (rank == 0) {
    // Now we have to create the output raster
    PRB_ERROR err = CreateOutputRaster(input_raster,
                                       conf.output_filename,
                                       input_raster->pixel_size(),
                                       conf.output_srs);

    if (err != PRB_NOERROR) {
      fprintf(stderr, "Error creating raster!\n");
    }
  }

  // Open the new output raster
  PTIFF* output_raster = open_raster(conf.output_filename);
  shared_ptr<ProjectedRaster> pr_output_raster(new ProjectedRaster(conf.output_filename));

  if (output_raster == NULL) {
    fprintf(stderr, "Could not open output raster\n");
    return 0;
  }

  // Now we will partition the output raster space. We will use a
  // maximum_height of 1 because we want single row or smaller
  // partitions to work with sptw.
  vector<Area> partitions = PartitionBySize(rank, 
					    process_count,
					    output_raster->y_size,
					    output_raster->x_size,
					    conf.partition_size,
					    1,
					    -1);

  // Now we loop through the returned partitions
  for (size_t i = 0; i < partitions.size(); ++i) {
    // Swap y axis
    Area swap = partitions.at(i);
    swap.ul.y = output_raster->y_size - partitions.at(i).ul.y - 1;
    swap.lr.y = output_raster->y_size - partitions.at(i).lr.y - 1;
    partitions.at(i) = swap;
    if (partitions.at(i).ul.x > partitions.at(i).lr.x) {
    printf("%f %f %f %f\n",partitions[i].ul.x, partitions[i].ul.y, partitions[i].lr.x,
	      partitions[i].lr.y);
    }
    // Calculate the input area used by the partition
    Area in_area = RasterMinbox(input_raster,
                                pr_output_raster,
                                partitions.at(i));

    // Read the input area 
    in_chunk = input_raster->create_raster_chunk(in_area);
    if (in_chunk == NULL) {
      fprintf(stderr, "Error reading input chunk! %f %f %f %f",
	      in_area.ul.x, in_area.ul.y, in_area.lr.x, in_area.lr.y);
      return 1;
    }

    // Create an allocated (not read from file) RasterChunk for the output 
    out_chunk = pr_output_raster->create_allocated_raster_chunk(partitions.at(i));
    if (out_chunk == NULL) {
      fprintf(stderr, "Error allocating output chunk! %f %f %f %f'n",
	      partitions[i].ul.x, partitions[i].ul.y, partitions[i].lr.x,
	      partitions[i].lr.y);
      return 1;
    }

    // Finally call ReprojectChunk
    bool ret = ReprojectChunk(in_chunk,
                              out_chunk,
                              "0",
                              "Min");

    delete in_chunk;
    delete out_chunk;
  }
              

  // Clean up
  close_raster(output_raster);
  output_raster = NULL;
  MPI_Finalize();

  return 0;
}
