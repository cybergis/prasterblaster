/*!
 * Copyright 0000 <Nobody>
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE
 *
 * This software is in the public domain, furnished "as is", without
 * technical support, and with no warranty, express or implied, as to
 * its usefulness for any purpose.
 *
 * @section DESCRIPTION
 *
 * Implementation of the Simple Parallel Tiff Writer
 *
 */

#include <endian.h>
#include <fcntl.h>
#include <gdal_priv.h>
#include <cpl_string.h>
#include <ogr_api.h>
#include <ogr_spatialref.h>
#include <mpi.h>
#include <sstream>
#include <tiff.h>
#include <tiffio.h>

#include "src/demos/sptw.h"
#include "src/rasterchunk.h"
#include "src/std_int.h"

using std::string;
using librasterblaster::RasterChunk;

namespace sptw {
SPTW_ERROR create_raster(string filename,
                         int64_t x_size,
                         int64_t y_size,
                         int band_count,
                         GDALDataType band_type,
                         double *geotransform,
                         string projection_srs) {
  GDALDriver *gtiff_driver = NULL;
  GDALDataset *ds = NULL;
  char **options = NULL;
  OGRSpatialReference srs;
  const char *format = "GTiff";

  GDALAllRegister();

  gtiff_driver = GetGDALDriverManager()->GetDriverByName(format);

  if (gtiff_driver == NULL) {
    return SP_CreateError;
  }

  options = CSLSetNameValue(options, "BIGTIFF", "YES");
  options = CSLSetNameValue(options, "INTERLEAVE", "PIXEL");
  options = CSLSetNameValue(options, "COMPRESS", "NONE");

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

  ds->SetGeoTransform(geotransform);

  // Close dataset
  GDALClose((GDALDatasetH) ds);

  return SP_None;
}

SPTW_ERROR create_tiled_raster(string filename,
                               int64_t x_size,
                               int64_t y_size,
                               int band_count,
                               GDALDataType band_type,
                               double *geotransform,
                               string projection_srs,
                               int64_t tile_size) {
  GDALDriver *gtiff_driver = NULL;
  GDALDataset *ds = NULL;
  char **options = NULL;

  GDALAllRegister();

  gtiff_driver = GetGDALDriverManager()->GetDriverByName("GTiff");

  if (gtiff_driver == NULL) {
    return SP_CreateError;
  }

  std::stringstream ts;
  ts << tile_size;

  options = CSLSetNameValue(options, "BIGTIFF", "YES");
  options = CSLSetNameValue(options, "INTERLEAVE", "PIXEL");
  options = CSLSetNameValue(options, "COMPRESS", "NONE");
  options = CSLSetNameValue(options, "TILED", "YES");
  options = CSLSetNameValue(options, "BLOCKXSIZE", ts.str().c_str());
  options = CSLSetNameValue(options, "BLOCKYSIZE", ts.str().c_str());

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

  ds->SetGeoTransform(geotransform);

  // Close dataset
  GDALClose((GDALDatasetH) ds);

  return SP_None;
}

PTIFF* open_raster(string filename) {
  PTIFF *ptiff = new PTIFF();
  char *c_filename = strdup(filename.c_str());

  GDALAllRegister();

  GDALDataset *ds = static_cast<GDALDataset*>(GDALOpen(filename.c_str(),
                                                       GA_Update));

  if (ds == NULL) {
    free(c_filename);
    return NULL;
  }

  ptiff->x_size = ds->GetRasterXSize();
  ptiff->y_size = ds->GetRasterYSize();
  ptiff->band_count = ds->GetRasterCount();
  ptiff->band_type = ds->GetRasterBand(1)->GetRasterDataType();
  ptiff->band_type_size = GDALGetDataTypeSize(ptiff->band_type)/8;
  ptiff->first_strip_offset = -1;
  ptiff->block_x_size = ptiff->x_size;
  ptiff->block_y_size = ptiff->y_size;

  GDALClose(ds);

  TIFF *tiffds = TIFFOpen(c_filename, "r");
  free(c_filename);

  if (tiffds == NULL) {
    fprintf(stderr, "Couldn't open tiff file\n");
    return NULL;
  }

  // Attempt to read TileWidth tag. If not found we assume file to have strips
  // instead.
  int64_t *offset = NULL;

  int ret = TIFFGetField(tiffds, TIFFTAG_TILEWIDTH, &(ptiff->block_x_size));
  ret = TIFFGetField(tiffds, TIFFTAG_TILELENGTH, &(ptiff->block_y_size));
  ret = TIFFGetField(tiffds, TIFFTAG_TILEOFFSETS, &offset);
  if (ret != 1) {
    ret = TIFFGetField(tiffds, TIFFTAG_STRIPOFFSETS, &offset);
    if (ret != 1) {
      fprintf(stderr, "Error reading strip offsets!\n");
      return NULL;
    } else {
      ptiff->first_strip_offset = *offset;
    }
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
    char *errstr = static_cast<char*>(malloc(5000));
    int errlen = 0;
    MPI_Error_string(rc, errstr, &errlen);
    fprintf(stderr,
            "MPI_File: Error opening file: %s: %s\n",
            c_filename,
            errstr);

    free(c_filename);
    free(errstr);
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
  delete ptiff;
  return SP_None;
}

SPTW_ERROR write_rows(PTIFF *ptiff,
                      void *buffer,
                      int64_t first_row,
                      int64_t last_row) {
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

SPTW_ERROR write_subrow(PTIFF *ptiff,
                        void *buffer,
                        int64_t row,
                        int64_t first_column,
                        int64_t last_column) {
  int64_t row_size = ptiff->x_size * ptiff->band_count * ptiff->band_type_size;
  int64_t subrow_size = (last_column - first_column + 1)
      * ptiff->band_count * ptiff->band_type_size;
  MPI_Offset offset = ptiff->first_strip_offset
      + (row * row_size) + first_column;

  MPI_Status status;

  MPI_File_write_at(ptiff->fh, offset, buffer, subrow_size, MPI_BYTE, &status);

  int count = 0;
  MPI_Get_count(&status, MPI_BYTE, &count);
  if (count != subrow_size) {
          fprintf(stderr,
                  "Error writing row! Wanted to write: %ld Actually wrote %d\n",
                  subrow_size, count);
  }

  return SP_None;
}

int64_t get_contiguous_size(PTIFF *tiff_file,
                            RasterChunk *chunk,
                            int64_t chunk_offset)
{
  const int64_t chunk_x_beginning = chunk->raster_location_.x + (chunk_offset % chunk->column_count_);
  const int64_t chunk_y_beginning = chunk->raster_location_.y + (chunk_offset / chunk->column_count_);
        
  const int64_t tile_x_beginning = (chunk_x_beginning / tiff_file->block_x_size) * tiff_file->block_x_size;
  const int64_t tile_y_beginning = (chunk_y_beginning / tiff_file->block_y_size) * tiff_file->block_y_size;
  const int64_t tile_x_end = tile_x_beginning + tiff_file->block_x_size - 1;
  const int64_t tile_y_end = tile_y_beginning + tiff_file->block_y_size - 1;
  const int64_t chunk_x_end = chunk->raster_location_.x + chunk->column_count_ - 1;
  const int64_t chunk_y_end = chunk->raster_location_.y + chunk->row_count_ - 1;

  librasterblaster::Coordinate end;
  end.x = tile_x_end;
  end.y = tile_y_end;

  if (chunk_x_end < tile_x_end) {
    //  The chunk is horizontally smaller than a tile so write a scanline.
    end.x = chunk_x_end;
    end.y = chunk_y_beginning;
  } else if (chunk_x_beginning != tile_x_beginning && chunk_x_end >= tile_x_end) {
    //  The chunk is wider than a tile so write the portion of the scanline that
    //  in in the tile.
    end.x = tile_x_end;
    end.y = chunk_y_beginning;
  } else if (chunk_x_beginning == tile_x_beginning
             && chunk_x_end == tile_x_end
             && tile_y_end > chunk_y_end) {
    //  The chunk is exactly the size of a tile so write the entire chunk.
    end.x = chunk_x_end;
    end.y = chunk_y_end;
  } else if (chunk_x_beginning == tile_x_beginning
             && chunk_x_end == tile_x_end
             && chunk_y_end >= tile_y_end) {
    //  The chunk is the width of a tile but vertically larger so write the rest
    //  of the tile.
    end.x = tile_x_end;
    end.y = tile_y_end;
  } else {
    return -1;
  }

  const int64_t first_chunk_x = chunk_offset % chunk->column_count_;
  const int64_t first_chunk_y = chunk_offset / chunk->column_count_;
  const int64_t last_chunk_x = end.x - chunk->raster_location_.x;
  const int64_t last_chunk_y = end.y - chunk->raster_location_.y;

  const int64_t contiguous_area_size = (last_chunk_x - first_chunk_x + 1) 
      + ((last_chunk_y - first_chunk_y) * chunk->column_count_);

  return contiguous_area_size;
}

int64_t chunk_to_file_offset(PTIFF *tiff_file,
                             RasterChunk *chunk,
                             int64_t chunk_offset) {
  //  First we have to translate the chunk_offset to raster coordinates.
  const int64_t raster_x = chunk->raster_location_.x 
      + (chunk_offset % chunk->column_count_);
  const int64_t raster_y = chunk->raster_location_.y
      + (chunk_offset / chunk->column_count_);
  const int64_t tile_padding = tiff_file->x_size % tiff_file->block_x_size != 0;
  const int64_t tile_columns = (tiff_file->x_size / tiff_file->block_x_size)
      + tile_padding;
  const int64_t tile_area = tiff_file->block_x_size 
      * tiff_file->block_y_size
      * tiff_file->band_count
      * tiff_file->band_type_size;
  const int64_t tile_row_size = tile_columns * tile_area;
  const int64_t tile_x_index = raster_x / tiff_file->block_x_size;
  const int64_t tile_y_index = raster_y / tiff_file->block_y_size;
  const int64_t tile_offset = (tile_x_index * tile_area)
      + tile_y_index * tile_row_size;


  const int64_t tile_x = raster_x % tiff_file->block_x_size;
  const int64_t tile_y = raster_y % tiff_file->block_y_size;


  return tile_offset + (tile_x * tiff_file->band_type_size)
      + (raster_y * tiff_file->block_x_size * tiff_file->band_type_size);
}

SPTW_ERROR write_rasterchunk(PTIFF *ptiff,
                             RasterChunk *chunk) {
  SPTW_ERROR err = SP_None;
  unsigned char *pixels = NULL;
  const int64_t x_offset = static_cast<int64_t>(chunk->raster_location_.x);
  const int64_t y_offset = static_cast<int64_t>(chunk->raster_location_.y);
  const int64_t chunk_size = chunk->row_count_
      * chunk->column_count_
      * chunk->band_count_
      * GDALGetDataTypeSize(chunk->pixel_type_)/8;
  int64_t chunk_offset = 0;
  int64_t file_offset = 0;
  int64_t write_size = 0;

  for (int64_t chunk_offset = 0; 
       chunk_offset < chunk_size;
       chunk_offset += write_size) {
    //  We need to calculate:
    //  1. offset into chunk->pixel buffer
    //  2. offset into file
    //  3. size of write from chunk->pixel buffer
    file_offset = ptiff->first_strip_offset + chunk_to_file_offset(ptiff,
                                                                   chunk,
                                                                   chunk_offset);
    write_size = get_contiguous_size(ptiff,
                                     chunk,
                                     chunk_offset);
    MPI_Status status;
    MPI_File_write_at(ptiff->fh,
                      file_offset,
                      &(static_cast<char*>(chunk->pixels_)[chunk_offset]),
                      write_size,
                      MPI_BYTE,
                      &status);
  }
  return SP_None;
}
}
