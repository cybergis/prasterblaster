/*!
 * Copyright 0000 <Nobody>
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
 * Header for high-level API
 *
 */
#ifndef SRC_DEMOS_PRASTERBLASTER_PIO_H_
#define SRC_DEMOS_PRASTERBLASTER_PIO_H_

#include "src/configuration.h"
#include "src/utils.h"

namespace librasterblaster {
PRB_ERROR prasterblasterpio(librasterblaster::Configuration conf);
}

#endif  // SRC_DEMOS_PRASTERBLASTER_PIO_H_

