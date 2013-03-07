///
/// Copyright 0000 <Nobody>
/// @file
/// @author David Matthew Mattli <dmattli@usgs.gov>
///
/// @section LICENSE
///
/// This software is in the public domain, furnished "as is", without
/// technical support, and with no warranty, express or implied, as to
/// its usefulness for any purpose.
///
/// @section DESCRIPTION
///
/// This file demonstrates how to use librasterblaster to implement parallel
/// raster reprojection. This implementation uses a MPI I/O to write to a tiff
/// file in parallel.
///
///

#include <vector>

#include "src/configuration.h"
#include "src/projectedraster.h"
#include "src/quadtree.h"
#include "src/reprojection_tools.h"
#include "src/sharedptr.h"

#include "src/demos/sptw.h"

using librasterblaster::Area;
using librasterblaster::PartitionBySize;
using librasterblaster::RasterChunk;
using librasterblaster::Configuration;
using librasterblaster::ProjectedRaster;
using librasterblaster::PRB_ERROR;
using librasterblaster::PRB_BADARG;
using librasterblaster::PRB_NOERROR;

using sptw::PTIFF;
using sptw::write_subrow;
using sptw::write_rasterchunk;
using sptw::open_raster;
using sptw::SPTW_ERROR;

void MyErrorHandler(CPLErr eErrClass, int err_no, const char *msg) {
        return;
}
/*! \page prasterblasterpio

\htmlonly
USAGE:
\endhtmlonly

\verbatim
prasterblasterpio [--t_srs target_srs] [--s_srs source_srs] 
                  [-r resampling_method] [-n partition_size]
                  [--dstnodata no_data_value]
                  source_file destination_file

\endverbatim

\section prasterblasterpio_description DESCRIPTION

<p> 

The prasterblasterpio demo program implements parallel raster reprojection and
demonstrates the use of librasterblaster. The implementation can be found in
prasterblaster-pio.cc.

</p>
 */

/** Main function for the prasterblasterpio program */
int prasterblaster_main(Configuration conf, int rank, int process_count) {
  RasterChunk *in_chunk, *out_chunk;

  // Replace CPLErrorHandler
  CPLPushErrorHandler(MyErrorHandler);

  GDALAllRegister();

  if (conf.input_filename == "" || conf.output_filename == "") {
    fprintf(stderr, "Specify an input and output filename\n");
    return 0;
  }

  // Open the input raster
  // 
  // The input raster is only read so we can use the serial i/o provided by the
  // GDAL library.
  GDALDataset *input_raster =
      static_cast<GDALDataset*>(GDALOpen(conf.input_filename.c_str(), GA_ReadOnly));
  if (input_raster == NULL) {
    fprintf(stderr, "Error opening input raster!\n");
    return 0;
  }

  // If we are the process with rank 0 we are responsible for the creation of
  // the output raster. 
  if (rank == 0) {
    printf("prasterblaster-pio: Beginning reprojection task\n");
    printf("\tInput File: %s, Output File: %s\n",
           conf.input_filename.c_str(), conf.output_filename.c_str());
    printf("\tInput Projection: %s\n\tOutput Projection: %s\n",
           input_raster->GetProjectionRef(), conf.output_srs.c_str());
    printf("\tProcess Count: %d\n", process_count);

    // Now we have to create the output raster
    printf("Creating output raster...\n");
    double gt[6];
    input_raster->GetGeoTransform(gt);
    PRB_ERROR err = librasterblaster::CreateOutputRaster(input_raster,
                                                         conf.output_filename,
                                                         gt[1],
                                                         conf.output_srs);
    if (err != PRB_NOERROR) {
      fprintf(stderr, "Error creating raster!: %d\n", err);
      return 1;
    }

    printf("Output raster created...\n");
  }

  // Wait for rank 0 to finish creating the file
  MPI_Barrier(MPI_COMM_WORLD);

  // Now open the new output file as a ProjectedRaster object. This object will
  // only be used to read metadata. It will _not_ be used to write to the output
  // file.
  PTIFF* output_raster = open_raster(conf.output_filename);
  GDALDataset *gdal_output_raster = 
      static_cast<GDALDataset*>(GDALOpen(conf.output_filename.c_str(), GA_ReadOnly));

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
                                            conf.partition_size);
  printf("Process %d of %d has %lu partitions with base size: %d\n",
         rank,
         process_count,
         static_cast<unsigned long>(partitions.size()),
         conf.partition_size);

  // Now we loop through the returned partitions
  for (size_t i = 0; i < partitions.size(); ++i) {
    // Swap y-axis of partition
    double t = partitions.at(i).lr.y;
    partitions.at(i).lr.y = partitions.at(i).ul.y;
    partitions.at(i).ul.y = t;

    // The RasterMinbox function calculates what part of the input raster
    // matches the given output partition.
    string input_wkt(input_raster->GetProjectionRef());
    string output_wkt(gdal_output_raster->GetProjectionRef());
    Area in_area = RasterMinbox(gdal_output_raster,
                                input_raster,
                                partitions.at(i));

    if (in_area.ul.x == -1.0) {  // chunk is outside of projected area
      // We should write the fill value here too
      continue;
    }

    // Now we use the ProjectedRaster object we created for the input file to
    // create a RasterChunk that has the pixel values read into it.
    in_chunk = RasterChunk::CreateRasterChunk(input_raster, in_area);

   PRB_ERROR chunk_err = RasterChunk::ReadRasterChunk(input_raster, in_chunk);
    if (chunk_err != PRB_NOERROR) {
      fprintf(stderr, "Error reading input chunk!\n");
      return 1;
    }

    // We want a RasterChunk for the output area but we area going to generate
    // the pixel values not read them from the file so we use
    // create_allocated_raster_chunk
    out_chunk = RasterChunk::CreateRasterChunk(gdal_output_raster,
                                               partitions.at(i));
    if (out_chunk == NULL) {
      fprintf(stderr, "Error allocating output chunk! %f %f %f %f'n",
              partitions[i].ul.x,
              partitions[i].ul.y,
              partitions[i].lr.x,
              partitions[i].lr.y);
      return 1;
    }
       
    // Now we call ReprojectChunk with the RasterChunk pair and the desired
    // resampler. ReprojectChunk performs the reprojection/resampling and fills
    // the output RasterChunk with the new values.
    bool ret = ReprojectChunk(in_chunk,
                              out_chunk,
                              conf.fillvalue,
                              conf.resampler);
    SPTW_ERROR err;
    err = write_rasterchunk(output_raster,
                            out_chunk);

    if (i % (partitions.size()/100) == 0) {
      printf("Rank: %d wrote chunk at %f %f, %d rows, %d columns\n\n", rank,
             out_chunk->raster_location_.x, out_chunk->raster_location_.y,
             out_chunk->row_count_, out_chunk->column_count_);
    }

    delete in_chunk;
    delete out_chunk;
  }

  // Clean up
  close_raster(output_raster);
  output_raster = NULL;

  return 0;
}

#if __USEMAIN

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

  int ret =  prasterblaster_main(conf, rank, process_count);
  MPI_Finalize();

  return ret;
}

#endif
