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
//
//

#include <gdal_priv.h>
#include <gdal.h>

#include "prasterblaster-pio.cc"

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "USAGE: rastercompare <control raster filename> "
            "<test raster filename>\n");
    return 1;
  }
  return librasterblaster::rastercompare(argv[1], argv[2]);
}

