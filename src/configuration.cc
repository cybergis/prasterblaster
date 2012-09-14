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

#include <getopt.h>

#include "src/configuration.h"

namespace librasterblaster {

struct option longopts[] = {
  {"output-projection", required_argument, NULL, 'p'},
  {"partition-count", required_argument, NULL, 'n'},
  {"resampler", required_argument, NULL, 'r'},
  {"fill-value", required_argument, NULL, 'f'},
  {"temporary-path", required_argument, NULL, 't'},
  {0, 0, 0, 0}
};

Configuration::Configuration() {
  resampler = NEAREST;
}

Configuration::Configuration(int argc, char *argv[]) {
  char c = 0;
  std::string arg = "";
  while ((c = getopt_long(argc, argv, "p:r:f:n:t:", longopts, NULL)) != -1) {
    switch (c) {
      case 0:
        // getopt_long() set a variable, just keep going
      case 'p':
        output_srs = optarg;
        break;
      case 'n':
        partition_count = strtol(optarg, NULL, 10);
        break;
      case 'r':
      arg = optarg;
        if (arg == "mean") {
          resampler = MEAN;
        } else if (arg == "nearest") {
          resampler = NEAREST;
        }
        break;
      case 'f':
        fillvalue = optarg;
        break;
      case 't':
        temporary_path = optarg;
        break;
      default:
        fprintf(stderr, "%s: option '-%c' is invalid: ignored\n",
                argv[0], optopt);
        break;
    }

    if (argc > 2) {
      output_filename = argv[argc-2];
      input_filename = argv[argc-1];
    }

    if (argc > 2 && temporary_path == "") {
      temporary_path = ".";
    }
  }

  return;
}
}
