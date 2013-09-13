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
#include "src/quadtree.h"
#include "src/reprojection_tools.h"

#include "src/demos/sptw.h"
#include "src/utils.h"

using librasterblaster::Area;
using librasterblaster::PartitionBySize;
using librasterblaster::RasterChunk;
using librasterblaster::Configuration;
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
namespace librasterblaster {

std::vector<Area> TilePartition(int rank,
                                int process_count,
                                PTIFF *tiff_file,
                                int tiles_per_partition) {
  std::vector<Area> parts = PartitionBySize(rank,
                                            process_count,
                                            tiff_file->tiles_down,
                                            tiff_file->tiles_across,
                                            tiles_per_partition);
  //  Now convert the tile partitions into raster coordinates
  for (int i = 0; i < parts.size(); ++i) {
    parts.at(i).ul.x *= tiff_file->block_x_size;
    parts.at(i).ul.y *= tiff_file->block_y_size;

    parts.at(i).lr.x += 1;
    parts.at(i).lr.x *= tiff_file->block_x_size;
    parts.at(i).lr.x -= 1;

    parts.at(i).lr.y += 1;
    parts.at(i).lr.y *= tiff_file->block_y_size;
    parts.at(i).lr.y -= 1;
  }

  return parts;
}

/** Main function for the prasterblasterpio program */
PRB_ERROR prasterblasterpio(Configuration conf) {
  RasterChunk *in_chunk, *out_chunk;
  double start_time, end_time;

  start_time = MPI_Wtime();

  int rank = 0;
  int process_count = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);

  // Replace CPLErrorHandler
  CPLPushErrorHandler(MyErrorHandler);

  GDALAllRegister();

  if (conf.input_filename == "" || conf.output_filename == "") {
    fprintf(stderr, "Specify an input and output filename\n");
    return PRB_BADARG;
  }

  // Open the input raster
  //
  // The input raster is only read so we can use the serial i/o provided by the
  // GDAL library.
  GDALDataset *input_raster =
      static_cast<GDALDataset*>(GDALOpen(conf.input_filename.c_str(),
                                         GA_ReadOnly));
  if (input_raster == NULL) {
    fprintf(stderr, "Error opening input raster!\n");
    return PRB_IOERROR;
  }

  // If we are the process with rank 0 we are responsible for the creation of
  // the output raster.
  if (rank == 0) {
    OGRSpatialReference sr;
    char *wkt;
    sr.SetFromUserInput(input_raster->GetProjectionRef());
    sr.exportToPrettyWkt(&wkt);
    printf("prasterblaster-pio: Beginning reprojection task\n");
    printf("\tInput File: %s, Output File: %s\n",
           conf.input_filename.c_str(), conf.output_filename.c_str());
    printf("\tInput Projection: %s\n\tOutput Projection: %s\n",
           wkt, conf.output_srs.c_str());
    printf("\tProcess Count: %d\n", process_count);
    OGRFree(wkt);

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
      return PRB_IOERROR;
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
      static_cast<GDALDataset*>(GDALOpen(conf.output_filename.c_str(),
                                         GA_ReadOnly));

  if (output_raster == NULL) {
    fprintf(stderr, "Could not open output raster\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    return PRB_IOERROR;
  }

  vector<Area> partitions;
  if (conf.partitioner == "tile") {
    partitions = TilePartition(rank,
                               process_count,
                               output_raster,
                               conf.partition_size);
  } else {
    partitions = PartitionBySize(rank,
                                 process_count,
                                 output_raster->y_size,
                                 output_raster->x_size,
                                 conf.partition_size);
  }
  if (rank == 0) {
    printf("Typical process has %lu partitions with base size: %d\n",
           static_cast<unsigned long>(partitions.size()),
           conf.partition_size);
  }
  double read_start, read_end, read_total;
  double write_start, write_end, write_total;
  double resample_start, resample_end, resample_total;
  double loop_start, prelude_end, prologue_end;

  read_total = write_total = resample_total = 0.0;

  // Now we loop through the returned partitions
  for (size_t i = 0; i < partitions.size(); ++i) {
    loop_start = MPI_Wtime();
    // Swap y-axis of partition
    double t = partitions.at(i).lr.y;
    partitions.at(i).lr.y = partitions.at(i).ul.y;
    partitions.at(i).ul.y = t;

    // Now we use the ProjectedRaster object we created for the input file to
    // create a RasterChunk that has the pixel values read into it.
    in_chunk = RasterChunk::CreateRasterChunk(input_raster,
                                              gdal_output_raster,
                                              partitions.at(i));
    prelude_end = MPI_Wtime();
    read_start = MPI_Wtime();
    PRB_ERROR chunk_err = RasterChunk::ReadRasterChunk(input_raster, in_chunk);
    if (chunk_err != PRB_NOERROR) {
      fprintf(stderr, "Error reading input chunk!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      return PRB_IOERROR;
    }
    read_end = MPI_Wtime();
    read_total += read_end - read_start;

    // We want a RasterChunk for the output area but we area going to generate
    // the pixel values not read them from the file so we use
    // CreateRasterChunk
    out_chunk = RasterChunk::CreateRasterChunk(gdal_output_raster,
                                               partitions.at(i));
    if (out_chunk == NULL) {
      fprintf(stderr, "Error allocating output chunk! %f %f %f %f\n",
              partitions[i].ul.x,
              partitions[i].ul.y,
              partitions[i].lr.x,
              partitions[i].lr.y);
      return PRB_BADARG;
    }

    // Now we call ReprojectChunk with the RasterChunk pair and the desired
    // resampler. ReprojectChunk performs the reprojection/resampling and fills
    // the output RasterChunk with the new values.
    resample_start = MPI_Wtime();
    bool ret = ReprojectChunk(in_chunk,
                              out_chunk,
                              conf.fillvalue,
                              conf.resampler);
    if (ret == false) {
            fprintf(stderr, "Error reprojecting chunk!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
    }
    resample_end = MPI_Wtime();
    resample_total += resample_end - resample_start;

    write_start = MPI_Wtime();
    SPTW_ERROR err;
    err = sptw::write_rasterchunk(output_raster,
                                  out_chunk);
    write_end = MPI_Wtime();
    write_total += write_end - write_start;
    if (err != sptw::SP_None) {
      fprintf(stderr, "Rank %d: Error writing chunk!\n", rank);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    delete in_chunk;
    delete out_chunk;
    prologue_end = MPI_Wtime();

    if (i % 10 == 0) {
      printf("<Rank %d wrote (%zd of %zd) in %f %f %f %f %f>  ",
             rank,
             i,
             partitions.size(),
             prelude_end - loop_start,
             read_end - read_start,
             write_end - write_start,
             resample_end - resample_start,
             prologue_end - loop_start);
      fflush(stdout);
    }
  }
  printf("Rank %d done\n", rank);
  // Clean up
  close_raster(output_raster);
  output_raster = NULL;
  delete gdal_output_raster;
  delete input_raster;

  // Report runtimes
  end_time = MPI_Wtime();
  double runtimes[4] = { end_time - start_time,
                         read_total,
                         write_total,
                         resample_total};
  std::vector<double> process_runtimes(process_count*4);
  MPI_Gather(runtimes,
             4,
             MPI_DOUBLE,
             &(process_runtimes[0]),
             4,
             MPI_DOUBLE,
             0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    printf("\nRank\t  Total Runtime\t  Read I/O\t  Write I/O\t  Resampling\n");
    for (unsigned int i = 0; i < process_runtimes.size(); i+=4) {
      printf("%4u\t| %.4fs\t| %.4fs\t| %.4fs\t| %.4fs\n",
             i/4,
             process_runtimes.at(i),
             process_runtimes.at(i+1),
             process_runtimes.at(i+2),
             process_runtimes.at(i+3));
    }
  }
  return PRB_NOERROR;
}
}
