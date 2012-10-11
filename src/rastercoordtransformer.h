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

/*
 *
 *
 */
#ifndef SRC_RASTERCOORDTRANSFORMER_H_
#define SRC_RASTERCOORDTRANSFORMER_H_

#include <ogr_spatialref.h>

#include "src/sharedptr.h"
#include "src/utils.h"



/// Raster Coordinate transformation class
/*
  This class implements the transformation of raster coordinates between two raster spaces with different projections and scales.
 */

namespace librasterblaster {
class ProjectedRaster;
class RasterCoordTransformer {
 public:
  // ! A constructor
  /* !  

    This constructor takes two initialized ProjectedRaster
    shared_ptrs and constructs a ready RasterCoordTransformer.
  */
  RasterCoordTransformer(shared_ptr<ProjectedRaster> source,
                         shared_ptr<ProjectedRaster> dest);

  // ! A constructor
  /* !

    This constructor takes a shared_ptr to a source raster and
    three parameters that describe a destination raster: a
    projection, upper-left coordinate, and a pixel size in
    meters.

  */
  RasterCoordTransformer(shared_ptr<ProjectedRaster> source,
                         shared_ptr<Projection> destination_projection,
                         Coordinate destination_ul,
                         double destination_pixel_size);

  // ! A constructor
  /* ! 

    This constructor takes six parameters. The first three
    describe the source raster: a projection, upper-left
    coordinate in projected coordinates, and a pixel size in meters.

    The second set of three parameters are the same as the first
    three but for the destination raster.

  */
  RasterCoordTransformer(shared_ptr<Projection> source_projection,
                         Coordinate source_ul,
                         double source_pixel_size,
                         int source_row_count,
                         int souce_column_count,
                         shared_ptr<Projection> destination_projection,
                         Coordinate destination_ul,
                         double destination_pixel_size);


  ~RasterCoordTransformer();

  // ! A normal member taking a single argument and returning an Area struct.
  /*

    This function takes a coordinate in the source raster space
    and maps is to an area in the destination raster space. The
    returned area consists of two points: an upper-left
    coordinate, and a lower-right. The coordinates are
    _inclusive_. So if the area consists of a single point, the
    upper-left == lower-right.

    \param source a Coordinate struct that specifies the point in the source raster space to map to the destination raster space.
  */
  Area Transform(Coordinate source, bool area_check = true);

 private:
  void init(shared_ptr<Projection> source_projection,
            Coordinate source_ul,
            double source_pixel_size,
            int source_row_count,
            int source_column_count,
            shared_ptr<Projection> destination_projection,
            Coordinate destination_ul,
            double destination_pixel_size);
  shared_ptr<Projection> src_proj, dest_proj;
  shared_ptr<OGRCoordinateTransformation> ctrans;
  Area maximum_geographic_area_;
  Coordinate source_ul_;
  double source_pixel_size_;
  Coordinate destination_ul_;
  double destination_pixel_size_;
};
}

#endif  // SRC_RASTERCOORDTRANSFORMER_H_
