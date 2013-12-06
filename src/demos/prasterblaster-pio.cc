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

#include <algorithm>
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
int simplerandom(int i) {
  return std::rand()%i;
}


int rastercompare(string control_filename, string test_filename) {
  const double delta = 0.001;
  GDALAllRegister();

  GDALDataset *control =
      static_cast<GDALDataset*>(GDALOpen(control_filename.c_str(), GA_ReadOnly));
  GDALDataset *test =
      static_cast<GDALDataset*>(GDALOpen(test_filename.c_str(), GA_ReadOnly));

  if (control == NULL) {
    fprintf(stderr, "Error opening control raster!\n");
    return 1;
  }

  if (test == NULL) {
    fprintf(stderr, "Error opening test raster!\n");
    return 1;
  }

  // Basic checks
  if (control->GetRasterYSize() != test->GetRasterYSize()) {
    printf("Control raster has %d rows and test has %d\n",
           control->GetRasterYSize(),
           test->GetRasterYSize());
    return 1;
  }

  if (control->GetRasterXSize() != test->GetRasterXSize()) {
    printf("Control raster has %d columns and test has %d\n",
           control->GetRasterXSize(),
           test->GetRasterXSize());
    return 1;
  }

  // Check pixel values
  const int band_count = control->GetRasterCount();
  GDALRasterBand *band = control->GetRasterBand(1);
  int block_x_size, block_y_size;
  band->GetBlockSize(&block_x_size, &block_y_size);

  double *control_pixels, *test_pixels;

  control_pixels = new double[band_count
                               * control->GetRasterXSize() * sizeof(double)];
  test_pixels = new double[band_count
                            * control->GetRasterXSize() * sizeof(double)];

  // Loop over rows of blocks
  //   Loop over columns of blocks
  //     Loop over rows in block
  //      Loop over columns in block
  int bad_pixels = 0;
  const int y_size = control->GetRasterYSize();
  const int x_size = control->GetRasterXSize();
  for (int y = 0; y < y_size; ++y) {
    // Read a row
    control->RasterIO(GF_Read,
                      0,
                      y,
                      control->GetRasterXSize(),
                      1,
                      static_cast<void*>(control_pixels),
                      control->GetRasterXSize(),
                      1,
                      GDT_Float64,
                      band_count,
                      NULL,
                      0,
                      0,
                      0);
    test->RasterIO(GF_Read,
                   0,
                   y,
                   test->GetRasterXSize(),
                   1,
                   static_cast<void*>(test_pixels),
                   test->GetRasterXSize(),
                   1,
                   GDT_Float64,
                   band_count,
                   NULL,
                   0,
                   0,
                   0);
    if (memcmp(control_pixels, test_pixels, x_size * sizeof(double)) != 0) {
      for (int x = 0; x < x_size; ++x) {
        if (fabs(control_pixels[x] - test_pixels[x]) > delta) {
          printf("Values at (%d, %d) are too different! %f vs %f\n",
                 x, y, control_pixels[x], test_pixels[x]);
          bad_pixels++;
        }
      }
    }
  }

  // Cleanup
  delete control_pixels;
  delete test_pixels;
  GDALClose(control);
  GDALClose(test);
  if (bad_pixels == 0) {
    return 0;
  } else {
    if (bad_pixels == 1) {
      printf("One pixel was different\n");
    } else {
      printf("%d pixels were different\n", bad_pixels);
    }
    return 1;
  }
}

std::vector<Area> TilePartition(int rank,
                                int process_count,
                                PTIFF *tiff_file,
                                int tiles_per_partition) {
  vector<Area> gparts;
  const int rows_per_partition = tiles_per_partition / tiff_file->tiles_across;
  int extra_rows = tiles_per_partition % tiff_file->tiles_across;
  const int partitions_per_row = tiff_file->tiles_across / tiles_per_partition;
  int extra_tiles = tiff_file->tiles_across % tiles_per_partition;
  printf("Beginning the loop!\n");
  if (rows_per_partition > 0) {
    printf("rows!\n");
    for (int i = 0; i < tiff_file->tiles_down; i += rows_per_partition) {
      int prows = rows_per_partition - 1;
      if (extra_rows > 0) {
        prows++;
        i++;
        extra_rows--;
      }
      gparts.push_back(Area(0, i, tiff_file->tiles_across-1, i+prows));
    }
  } else {
    printf("columns!\n");
    for (int i = 0; i < tiff_file->tiles_down; ++i) {
      for (int j = 0, tilei = 0; j < process_count; ++j) {
        int ptiles = tiles_per_partition - 1;
        if (j / extra_tiles == 0) {
          ptiles += extra_tiles % process_count;
        }
        gparts.push_back(Area(tilei, i, tilei+ptiles, i));
        tilei += ptiles + 1;
      }
    }

  }
  printf("BEAT THE LOOP!\n\n");
  std::vector<Area> partitions;
  // Sort the partition for better file locality
  std::sort(gparts.begin(), gparts.end(), partition_compare);

  // Now we shuffle the partitions for load balancing.
  // All processes should generate the same shuffle.
  for (size_t i = 0; i < gparts.size()-process_count; i += process_count) {
    vector<Area>::iterator beg = gparts.begin() + i;
    vector<Area>::iterator end = beg + process_count - 1;
    std::random_shuffle(beg, end, simplerandom);
  }

  for (size_t i = 0; i < gparts.size(); ++i) {
    if (i % process_count == rank) {
      partitions.push_back(gparts.at(i));
    }
  }

  //  Now convert the tile partitions into raster coordinates
  for (int i = 0; i < partitions.size(); ++i) {
    partitions.at(i).ul.x *= tiff_file->block_x_size;
    partitions.at(i).ul.y *= tiff_file->block_y_size;

    partitions.at(i).lr.x += 1;
    partitions.at(i).lr.x *= tiff_file->block_x_size;
    partitions.at(i).lr.x -= 1;

    partitions.at(i).lr.y += 1;
    partitions.at(i).lr.y *= tiff_file->block_y_size;
    partitions.at(i).lr.y -= 1;
  }

  return partitions;
}

/** Main function for the prasterblasterpio program */
PRB_ERROR prasterblasterpio(Configuration conf) {
  RasterChunk *in_chunk, *out_chunk;
  double start_time, end_time, preloop_time;

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
  double read_total, misc_start, misc_total;
  double write_start, write_end, write_total;
  double resample_start, resample_end, resample_total;
  double loop_start, prelude_end, minbox_total;

  read_total = write_total = resample_total = misc_total = minbox_total = 0.0;
  preloop_time = MPI_Wtime() - start_time;

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
    minbox_total += MPI_Wtime() - loop_start;

    PRB_ERROR chunk_err = RasterChunk::ReadRasterChunk(input_raster, in_chunk);
    if (chunk_err != PRB_NOERROR) {
      fprintf(stderr, "Error reading input chunk!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      return PRB_IOERROR;
    }
    read_total += MPI_Wtime() - prelude_end;

    misc_start = MPI_Wtime();
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
    misc_total += MPI_Wtime() - misc_start;

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
    if (err != sptw::SP_None) {
      fprintf(stderr, "Rank %d: Error writing chunk!\n", rank);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    write_end = MPI_Wtime();
    write_total += write_end - write_start;

    misc_start = MPI_Wtime();
    delete in_chunk;
    delete out_chunk;

    if (i % 10 == 0) {
      printf("<Rank %d wrote (%zd of %zd)> ",
             rank,
             i,
             partitions.size());
      fflush(stdout);
    }
    misc_total += MPI_Wtime() - misc_start;
  }
  printf("Rank %d done\n", rank);
  // Clean up
  write_start = MPI_Wtime();
  close_raster(output_raster);
  write_total += MPI_Wtime() - write_start;

  misc_start = MPI_Wtime();
  output_raster = NULL;
  delete gdal_output_raster;
  delete input_raster;
  misc_total += MPI_Wtime() - misc_start;

  // Report runtimes
  end_time = MPI_Wtime();
  double runtimes[7] = { end_time - start_time,
                         preloop_time,
                         minbox_total,
                         read_total,
                         resample_total,
                         write_total,
                         misc_total };
  std::vector<double> process_runtimes(process_count*7);
  MPI_Gather(runtimes,
             7,
             MPI_DOUBLE,
             &(process_runtimes[0]),
             7,
             MPI_DOUBLE,
             0,
             MPI_COMM_WORLD);

  FILE *timing_file = stdout;
  if (rank == 0) {
    if (conf.timing_filename != "") {
      timing_file = fopen(conf.timing_filename.c_str(), "w");
      if (timing_file == NULL) {
        fprintf(stderr, "Error creating timing output file");
        timing_file = stdout;
      }
    }

    fprintf(timing_file, "process,total,preloop,minbox,read,resample,write,misc\n");
    for (unsigned int i = 0; i < process_runtimes.size(); i+=7) {
      fprintf(timing_file, "%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
              i/7,
              process_runtimes.at(i),
              process_runtimes.at(i+1),
              process_runtimes.at(i+2),
              process_runtimes.at(i+3),
              process_runtimes.at(i+4),
              process_runtimes.at(i+5),
              process_runtimes.at(i+6));
    }
  }

  if (rank == 0 && conf.timing_filename != "") {
    fclose(timing_file);
  }

  return PRB_NOERROR;
}
}
