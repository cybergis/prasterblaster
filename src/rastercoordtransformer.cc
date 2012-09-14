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
//
/* ! \mainpage RasterCoordTransformer
 *
 *
 */

#include <ogr_api.h>
#include <ogr_spatialref.h>
#include <gdal.h>
#include <gdal_priv.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/coordinate.h"

#include "projectedraster.h"
#include "reprojection_tools.h"
#include "resampler.h"
#include "sharedptr.h"
#include "quadtree.h"

namespace librasterblaster {

RasterCoordTransformer::
RasterCoordTransformer(shared_ptr<ProjectedRaster> source,
                       shared_ptr<ProjectedRaster> _dest) {
  src_proj = shared_ptr<Projection>(source->projection());
  source_pixel_size_ = source->pixel_size();
  source_ul_ = Coordinate(source->ul_x(), source->ul_y(), UNDEF);
  dest_proj = shared_ptr<Projection>(_dest->projection());
  destination_pixel_size_ = _dest->pixel_size();
  destination_ul_ = Coordinate(_dest->ul_x(), _dest->ul_y(), UNDEF);
  return;
}

RasterCoordTransformer::
RasterCoordTransformer(shared_ptr<ProjectedRaster> source,
                       shared_ptr<Projection> destination_projection,
                       Coordinate destination_ul,
                       double destination_pixel_size) {
  destination_pixel_size_ = destination_pixel_size;
  destination_ul_ = destination_ul;
  source_ul_ = Coordinate(source->ul_x(), source->ul_y(), UNDEF);
  source_pixel_size_ = source->pixel_size();

  dest_proj = destination_projection;
  src_proj =  shared_ptr<Projection>(source->projection());

  return;
}

RasterCoordTransformer::
RasterCoordTransformer(shared_ptr<Projection> source_projection,
                       Coordinate source_ul,
                       double source_pixel_size,
                       shared_ptr<Projection> destination_projection,
                       Coordinate destination_ul,
                       double destination_pixel_size) {
  src_proj = source_projection;
  source_ul_ = source_ul;
  source_pixel_size_ = source_pixel_size;
  dest_proj = destination_projection;
  destination_ul_ = destination_ul;
  destination_pixel_size_ = destination_pixel_size;
  return;
}

RasterCoordTransformer::~RasterCoordTransformer() {
  return;
}

Area RasterCoordTransformer::
Transform(Coordinate source) {
  Area value;
  Coordinate temp1, temp2;

  temp1.x = temp1.y = 0.0;

  value.ul = temp1;
  value.lr = temp1;

  temp1.x = (static_cast<double>(source.x) * source_pixel_size_) + source_ul_.x;
  temp1.y = (static_cast<double>(source.y) * source_pixel_size_) - source_ul_.y;

  src_proj->inverse(temp1.x, temp1.y, &temp2.x, &temp2.y);
  src_proj->forward(temp2.x, temp2.y, &temp2.x, &temp2.y);

  if (fabs(temp1.x - temp2.x) > 0.01) {
    // Point is outside defined projection area, return no-value
    value.ul.x = -1.0;
    value.lr.x = -1.0;
    return value;
  }

  temp1.x = (static_cast<double>(source.x) * source_pixel_size_) + source_ul_.x;
  temp1.y = (static_cast<double>(source.y) * source_pixel_size_) - source_ul_.y;
  temp2 = temp1;

  // Now we are going to assign temp1 as the UL of our pixel and
  // temp2 as LR
  temp2.x += source_pixel_size_/2;
  temp2.y -= source_pixel_size_/2;

  src_proj->inverse(temp1.x, temp1.y, &temp1.x, &temp1.y);
  dest_proj-> forward(temp1.x, temp1.y, &temp1.x, &temp1.y);
  src_proj->inverse(temp2.x, temp2.y, &temp2.x, &temp2.y);
  dest_proj-> forward(temp2.x, temp2.y, &temp2.x, &temp2.y);
  // temp1/temp2 now contain coords to input projection
  // Now convert to points in the raster coordinate space.
  temp1.x -= destination_ul_.x;
  temp1.y += destination_ul_.y;
  temp1.x /= destination_pixel_size_;
  temp1.y /= destination_pixel_size_;
  temp2.x -= destination_ul_.x;
  temp2.y += destination_ul_.y;
  temp2.x /= destination_pixel_size_;
  temp2.y /= destination_pixel_size_;

  value.ul = temp1;
  value.lr = temp2;

  return value;
}
}
