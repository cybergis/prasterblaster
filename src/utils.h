/**
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
 * Helper utilities and definitions
 *
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include "src/gctp_cpp/coordinate.h"

namespace librasterblaster {
/** An enum of possible return values */
enum PRB_ERROR {
  PRB_NOERROR, /*!< No error occurred */
  PRB_IOERROR, /*!< Error communicating or performing file I/O */
  PRB_BADARG,  /*!< Bad argument provided */
  PRB_PROJERROR, /*!< Error with projection specification */
};

/** 
 * @class 
 */
struct Area {
  /// A constructor
  /** 
   * @brief This constructor initializes points to zero.
   */
  Area() {}
  /// A constructor
  /**
   * @brief This constructor allows both coordinates to be given initial values.
   * @param ulx Value for the upper-left x 
   * @param uly Value for the upper-left y
   * @param ulx Value for the lower-right x 
   * @param ulx Value for the lower-right y
   */
Area(double ulx,
     double uly,
     double lrx,
     double lry) : ul(ulx, uly, UNDEF), lr(lrx, lry, UNDEF) {}
  /// Upper-left coordinate of area
  Coordinate ul;
  /// Lower-right coordinate of area
  Coordinate lr;
  /// Units used in ul and lr coordinates
  ProjUnit units;
};
}

#endif  // SRC_UTILS_H_
