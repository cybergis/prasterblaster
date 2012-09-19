
#include <endian.h>
#include <fcntl.h>
#include <gdal_priv.h>
#include <cpl_string.h>
#include <ogr_api.h>
#include <ogr_spatialref.h>
#include <mpi.h>
#include <stdint.h>
#include <tiff.h>
#include <tiffio.h>

#include "sptw.h"

using std::string;

namespace sptw {
SPTW_ERROR create_raster(string filename,
                         int x_size,
                         int y_size,
                         double ul_x,
                         double ul_y,
                         double pixel_size,
                         int band_count,
                         GDALDataType band_type,
                         string projection_srs) {
  GDALDriver *gtiff_driver = NULL;
  GDALDataset *ds = NULL;
  char **options = NULL;
  OGRSpatialReference srs;
  char *srs_wkt = NULL;
  GDALRasterBand *band = NULL;
  const char *format = "GTiff";
  double geotransform[16] = {0.0};

  GDALAllRegister();

  gtiff_driver = GetGDALDriverManager()->GetDriverByName(format);

  if (gtiff_driver == NULL) {
    return SP_CreateError;
  }  

  options = CSLSetNameValue( options, "BIGTIFF", "YES" );

  ds = gtiff_driver->Create(filename.c_str(),
                            x_size,
                            y_size,
                            band_count,
                            band_type,
                            options);

  // Clean up options
  CSLDestroy(options);

  CPLErr err = ds->SetProjection(projection_srs.c_str());

  if (err != CE_None) {
    return SP_BadArg;
  }

  geotransform[1] = geotransform[5] = pixel_size;
  geotransform[0] = ul_x;
  geotransform[3] = ul_y;

  ds->SetGeoTransform(geotransform);
  

  // Close dataset
  GDALClose((GDALDatasetH) ds);

  return SP_None;
}

PTIFF* open_raster(string filename) {
  PTIFF *ptiff = new PTIFF();
  char *c_filename = strdup(filename.c_str());
  
  GDALAllRegister();
  
  GDALDataset *ds = (GDALDataset*)GDALOpen(filename.c_str(), GA_Update);

  if (ds == NULL) {
    return NULL;
  }
  
  ptiff->x_size = ds->GetRasterXSize();
  ptiff->y_size = ds->GetRasterYSize();
  ptiff->band_count = ds->GetRasterCount();
  ptiff->band_type = ds->GetRasterBand(1)->GetRasterDataType();
  ptiff->band_type_size = GDALGetDataTypeSize(ptiff->band_type)/8;

  GDALClose(ds);

  TIFF *tiffds = TIFFOpen(c_filename, "r");
  free(c_filename);

  if (tiffds == NULL) {
    fprintf(stderr, "Couldn't open tiff file\n");
    return NULL;
  }

  uint64_t *offset = NULL;

  int ret = TIFFGetField(tiffds,TIFFTAG_STRIPOFFSETS, &offset);
  if (ret != 1) {
    fprintf(stderr, "Error reading strip offsets!\n");
  }

  if (offset == NULL) {
    fprintf(stderr, "Error reading strip offsets\n");
    free(c_filename);
    return NULL;
  }
  
  ptiff->first_strip_offset = *offset;

  TIFFClose(tiffds);

  c_filename = strdup(filename.c_str());
  int rc = MPI_File_open(MPI_COMM_WORLD,
			 c_filename,
                         MPI_MODE_RDWR,
                         MPI_INFO_NULL,
                         &(ptiff->fh));

  if (rc != MPI_SUCCESS) {
    char *errstr = (char*)malloc(5000);
    int errlen = 0;
    MPI_Error_string(rc, errstr, &errlen);
    fprintf(stderr, "MPI_File: Error opening file: %s: %s\n", c_filename, errstr);

    return NULL;
  }
  
  MPI_File_set_atomicity(ptiff->fh, 0);

  if (c_filename != NULL) {
    free(c_filename);
  }
  
  return ptiff;
}

SPTW_ERROR close_raster(PTIFF *ptiff) {
  MPI_File_close(&(ptiff->fh));

  return SP_None;
} 

SPTW_ERROR write_rows(PTIFF *ptiff, char *buffer, int64_t first_row, int64_t last_row) {

  int64_t row_size = ptiff->x_size * ptiff->band_count * ptiff->band_type_size;
  MPI_Offset offset = ptiff->first_strip_offset + (first_row * row_size);
  int count = (last_row - first_row + 1) * row_size;

  MPI_Status status;
  int err = MPI_SUCCESS;

  err = MPI_File_write_at(ptiff->fh, offset, buffer, count, MPI_BYTE, &status);

  if (err != MPI_SUCCESS) {
    fprintf(stderr, "Write error: %d\n", status.MPI_ERROR);
    return SP_WriteError;
  }

  return SP_None;
}

SPTW_ERROR write_subrow(PTIFF *ptiff, char *buffer, int64_t row, int64_t first_column, int64_t last_column) {
  int64_t row_size = ptiff->x_size * ptiff->band_count * ptiff->band_type_size;
  int64_t subrow_size = (first_column - last_column + 1) * ptiff->band_count * ptiff->band_type_size;
  MPI_Offset offset = ptiff->first_strip_offset + (row * row_size);

  MPI_Status status;

  MPI_File_write_at(ptiff->fh, offset, buffer, subrow_size, MPI_BYTE, &status);

  if (status.MPI_ERROR != MPI_SUCCESS) {
    return SP_WriteError;
  }

  return SP_None;
}
}