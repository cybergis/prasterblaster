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
// The RasterChunk class represents an in-memory, georeferenced section of
// raster.
//
//


#ifndef SRC_RASTERCHUNK_H_
#define SRC_RASTERCHUNK_H_

#include <string>

#include <gdal_priv.h>

#include "src/gctp_cpp/projection.h"
#include "src/gctp_cpp/coordinate.h"
#include "src/sharedptr.h"
#include "src/utils.h"

namespace librasterblaster {
/// A class representing an in-memory part of a raster.
class RasterChunk {
 public:
  static RasterChunk* CreateRasterChunk(GDALDataset *ds, Area chunk_area);
  static PRB_ERROR ReadRasterChunk(GDALDataset *ds, RasterChunk *chunk);
  static PRB_ERROR WriteRasterChunk(GDALDataset *ds, RasterChunk *chunk);
  /// RasterChunk constructor
  RasterChunk() {
    this->pixels_ = NULL;
  }
  /// RasterChunk destructor
  /**
   * This destructor frees the memory, if any,  allocated at pixels_.
   */
  ~RasterChunk() {
    if (this->pixels_ != NULL) {
      free(this->pixels_);
    }
  }
  std::string projection_;
  /// Location of the chunk, in raster coordinates
  /** 
   * This variable represents the upper-left location of the raster chunk, in
   * the raster coordinates.
   */
  Coordinate raster_location_;
  /// Upper-left corner, in projected coordinates
  /**
   * This variable represents the upper-left location of the raster chunk in
   * projected coordinates.
   */
  Coordinate ul_projected_corner_;
  /// Size of pixel, in meters
  /**
   * This variable represents the size of the pixel in meters.
   */
  double pixel_size_;  // in meters
  /// Number of rows
  int row_count_;
  /// Number of columns
  int column_count_;
  /// Datatype of pixel values
  GDALDataType pixel_type_;
  /// Number of bands
  int band_count_;
  /// GDAL geotransform
  double geotransform_[6];
  /// Pointer to pixel values
  void *pixels_;
};
}


#endif  // SRC_RASTERCHUNK_H_
