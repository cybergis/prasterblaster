/*!
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE
 * This work was produced as a part of the official duties of a
 * federal employee and is in the public domain.

 * @section DESCRIPTION
 *
 *
 *
 */


#ifndef RASTERCHUNK_HH_
#define RASTERCHUNK_HH_


#include <gdal_priv.h>


#include "projectedraster.hh"

namespace RasterChunk {

class ChunkExtent
{
  int first_index;
  int last_index;
};

class RasterChunk 
{
  ChunkExtent extent;
  GDALDataType *pixel_type;
  int band_count;
  double geotransform[6];
  void *pixels;
};

}

#endif //RASTERCHUNK_HH_
