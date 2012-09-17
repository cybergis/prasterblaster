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
// The Configuration class represents the configuration of a reprojection task.
//
//


#ifndef SRC_LIBRASTERBLASTER_CONFIGURATION_H_
#define SRC_LIBRASTERBLASTER_CONFIGURATION_H_

#include <string>

#include "src/resampler.h"
#include "src/reprojection_tools.h"

using std::string;

namespace librasterblaster {
struct Configuration {
  Configuration();
  Configuration(int argc, char *argv[]);

  string input_filename;
  string input_srs;
  string output_filename;
  string output_srs;

  RESAMPLER resampler;
  string fillvalue;
  string nodata_value;
  string temporary_path;
  long partition_count;

};
}

#endif  // SRC_LIBRASTERBLASTER_CONFIGURATION_H
