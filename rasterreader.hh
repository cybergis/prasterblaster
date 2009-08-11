
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
