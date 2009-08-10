
#ifndef RASTER_READER_HH
#define RASTER_READER_HH

#include <string>

#include "projectedraster.hh"

class RasterReader
{
public:
  static ProjectedRaster* readRaster(std::string filename);
  // filename is without extensionx
  static ProjectedRaster* readImgRaster(std::string filename);
  static void writeRaster(std::string filename,
			  ProjectedRaster *raster);

};


#endif // RASTER_READER_HH
