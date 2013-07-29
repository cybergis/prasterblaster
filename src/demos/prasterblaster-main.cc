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
// Implements main() function for prasterblasterpio command-line tool.
//
//

#include <mpi.h>

#include "src/configuration.h"
#include "src/demos/prasterblaster-pio.h"

using librasterblaster::Configuration;
using librasterblaster::prasterblasterpio;
using librasterblaster::RasterChunk;

int main(int argc, char *argv[]) {
  int rank = 0;
  int process_count = 0;
  RasterChunk *in_chunk, *out_chunk;

  // Give MPI_Init first run at the command-line arguments
  MPI_Init(&argc, &argv);

  // Initialize Configuration object
  Configuration conf(argc, argv);

  int ret =  prasterblasterpio(conf);
  MPI_Finalize();

  return ret;
}
