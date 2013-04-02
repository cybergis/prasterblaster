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

#include <cmath>

#include "src/gctp_cpp/projection.h"
#include "src/gctp_cpp/transformer.h"
#include "src/gctp_cpp/coordinate.h"

#include "src/projectedraster.h"
#include "src/reprojection_tools.h"
#include "src/resampler.h"
#include "src/sharedptr.h"
#include "src/quadtree.h"

namespace librasterblaster {
shared_ptr<Projection> srs_to_projection(string projection_string) {
  OGRSpatialReference srs;
  shared_ptr<Projection> out_proj;
  char *wkt = NULL;
  char *tmp;

  OGRErr err = srs.importFromProj4(projection_string.c_str());
  if (err != OGRERR_NONE) {
    wkt = strdup(projection_string.c_str());
    tmp = wkt;
    err = srs.importFromWkt(&tmp);
    free(wkt);
  }

  if (err != OGRERR_NONE) {
    fprintf(stderr, "Error parsing projection: %s\n",
        projection_string.c_str());
    return out_proj;
  }

  int64_t proj_code, datum_code, zone;
  double *params = NULL;

  srs.exportToUSGS(&proj_code, &zone, &params, &datum_code);

  out_proj = shared_ptr<Projection>(
      Transformer::convertProjection(static_cast<ProjCode>(proj_code)));

  if (!out_proj) {
    return out_proj;
  }

  out_proj->setUnits(METER);
  out_proj->setDatum(static_cast<ProjDatum>(datum_code));
  out_proj->setParams(params);

  OGRFree(params);

  return out_proj;
}

RasterCoordTransformer::
RasterCoordTransformer(string source_projection,
                       Coordinate source_ul,
                       double source_pixel_size,
                       int source_row_count,
                       int source_column_count,
                       string destination_projection,
                       Coordinate destination_ul,
                       double destination_pixel_size) {
  init(source_projection,
       source_ul,
       source_pixel_size,
       source_row_count,
       source_column_count,
       destination_projection,
       destination_ul,
       destination_pixel_size);
  return;
}

RasterCoordTransformer::~RasterCoordTransformer() {
  return;
}

void RasterCoordTransformer::init(string source_projection,
                                  Coordinate source_ul,
                                  double source_pixel_size,
                                  int source_row_count,
                                  int source_column_count,
                                  string destination_projection,
                                  Coordinate destination_ul,
                                  double destination_pixel_size) {
  src_proj = srs_to_projection(source_projection);
  source_ul_ = source_ul;
  source_pixel_size_ = source_pixel_size;
  dest_proj = srs_to_projection(destination_projection);
  destination_ul_ = destination_ul;
  destination_pixel_size_ = destination_pixel_size;

  OGRSpatialReference source_sr, dest_sr;
  char *source_wkt = strdup(source_projection.c_str());
  char *dest_wkt = strdup(destination_projection.c_str());

  source_sr.SetFromUserInput(source_projection.c_str());
  dest_sr.SetFromUserInput(destination_projection.c_str());

  OGRCoordinateTransformation *t = OGRCreateCoordinateTransformation(&source_sr,
                                                                     &dest_sr);
  free(source_wkt);
  free(dest_wkt);

  if (t != NULL) {
    ctrans.reset(t);
  } else {
    printf("BAD!\n\n");
    return;
  }

  Coordinate ul, lr;
  ul = source_ul_;
  src_proj->inverse(ul.x, ul.y, &ul.x, &ul.y);
  lr.x = source_ul_.x + (source_pixel_size_ * source_column_count);
  lr.y = source_ul_.y - (source_pixel_size_ * source_row_count);
  src_proj->inverse(lr.x, lr.y, &lr.x, &lr.y);
  ul.x = -179.99;
  lr.x = 179.99;
  maximum_geographic_area_.ul = ul;
  maximum_geographic_area_.lr = lr;
  return;
}

Area RasterCoordTransformer::
Transform(Coordinate source, bool area_check) {
  Area value;
  Coordinate temp1, temp2;

  temp1.x = temp1.y = 0.0;

  value.ul = temp1;
  value.lr = temp1;

  temp1.x = (static_cast<double>(source.x) * source_pixel_size_) + source_ul_.x;
  temp1.y = source_ul_.y - (static_cast<double>(source.y) * source_pixel_size_);

  src_proj->inverse(temp1.x, temp1.y, &temp2.x, &temp2.y);
  src_proj->forward(temp2.x, temp2.y, &temp2.x, &temp2.y);

  if (area_check && fabs(temp1.x - temp2.x) > 0.01) {
    // Point is outside defined projection area, return no-value
    value.ul.x = -1.0;
    value.lr.x = -1.0;
    return value;
  }

  temp1.x = (static_cast<double>(source.x) * source_pixel_size_) + source_ul_.x;
  temp1.y = source_ul_.y - (static_cast<double>(source.y) * source_pixel_size_);
  temp2 = temp1;

  // Now we are going to assign temp1 as the UL of our pixel and
  // temp2 as LR
  temp2.x += sqrt(2 * source_pixel_size_ * source_pixel_size_);
  temp2.y -= sqrt(2 * source_pixel_size_ * source_pixel_size_);

//  src_proj->inverse(temp1.x, temp1.y, &temp1.x, &temp1.y);
//  dest_proj->forward(temp1.x, temp1.y, &temp1.x, &temp1.y);

  ctrans->Transform(1, &temp1.x, &temp1.y);
//  src_proj->inverse(temp2.x, temp2.y, &temp2.x, &temp2.y);
//  dest_proj->forward(temp2.x, temp2.y, &temp2.x, &temp2.y);
  ctrans->Transform(1, &temp2.x, &temp2.y);

  // temp1/temp2 now contain coords to input projection
  // Now convert to points in the raster coordinate space.
  temp1.x -= destination_ul_.x;
  temp1.y = destination_ul_.y - temp1.y;
  temp1.x /= destination_pixel_size_;
  temp1.y /= destination_pixel_size_;
  temp2.x -= destination_ul_.x;
  temp2.y = destination_ul_.y - temp2.y;
  temp2.x /= destination_pixel_size_;
  temp2.y /= destination_pixel_size_;

  value.ul = temp1;
  value.lr = temp2;

  // Now validate and round pixel values
  // Truncate values
  value.ul.x = floor(fabs(value.ul.x));
  value.ul.y = floor(fabs(value.ul.y));
  value.lr.x = floor(fabs(value.lr.x));
  value.lr.y = floor(fabs(value.lr.y));

  // Check that entries are valid
  if (value.ul.x < 0.0 ||
      value.lr.x < 0.0) {
    value.ul.x = -1.0;
    value.lr.x = -1.0;
  }

  if (value.ul.x > value.lr.x) {
    value.lr.x = value.ul.x;
  }

  if (value.ul.y > value.lr.y) {
    value.lr.y = value.ul.y;
  }

  return value;
}

bool RasterCoordTransformer::ready() {
  if (ctrans.get() == NULL) {
    return false;
  }
  return true;
}
}
