


#ifndef SRC_SPTW_H_
#define SRC_SPTW_H_

#include <gdal_priv.h>
#include <mpi.h>
#include <stdint.h>

using std::string;

namespace sptw {
enum SPTW_ERROR {
  SP_None,
  SP_CreateError,
  SP_WriteError,
  SP_BadArg,
};

struct PTIFF {
  MPI_File fh;
  int x_size;
  int y_size;
  int band_count;
  GDALDataType band_type;
  int band_type_size;
  uint32_t first_strip_offset;
};

SPTW_ERROR create_raster(string filename,
                         int x_size,
                         int y_size,
                         int band_count,
                         double ul_x,
                         double ul_y,
                         double pixel_size,
                         GDALDataType band_type,
                         string projection_srs);

PTIFF* open_raster(string filename);
SPTW_ERROR close_raster(PTIFF *ptiff);
SPTW_ERROR write_rows(PTIFF *ptiff, 
		      void *buffer, 
		      int64_t first_row, 
		      int64_t last_row);
SPTW_ERROR write_subrow(PTIFF *ptiff, 
			void *buffer, 
			int64_t row, 
			int64_t first_column, 
			int64_t last_column);
}

#endif // SRC_SPTW_H_
