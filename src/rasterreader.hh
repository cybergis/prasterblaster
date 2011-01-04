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
 */


#ifndef RASTER_READER_HH
#define RASTER_READER_HH

#include <string>

#include "projectedraster.hh"

class RasterReader
{
public:
  static ProjectedRaster* readRaster(std::string filename);
  // filename is without extensionx
  static ProjectedRaster* readImgRaster(std::string filename, int rank = 0, int num_procs = 1);
  static void writeRaster(std::string filename,
			  ProjectedRaster *raster);
  static ProjectedRaster* readSubraster(std::string filename, int rank, int num_procs);

};


#endif // RASTER_READER_HH
