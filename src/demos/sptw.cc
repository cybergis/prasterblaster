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

#include <algorithm>
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
#include <vector>

#include "src/demos/sptw.h"
#include "src/rasterchunk.h"
#include "src/std_int.h"
#include "src/utils.h"

using std::string;
using librasterblaster::RasterChunk;
using librasterblaster::Area;

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


  return tiff_file->first_strip_offset + tile_offset + (tile_x * tiff_file->band_type_size)
      + (tile_y * tiff_file->block_x_size * tiff_file->band_type_size);
}

SPTW_ERROR fill_stack(std::vector<Area> *write_stack,
                      Area old_area,
                      Area written_subset) {
  const double size_above = written_subset.ul.y - old_area.ul.y;
  const double size_below = old_area.lr.y - written_subset.lr.y;
  const double size_left = written_subset.ul.x - old_area.ul.x;
  const double size_right = old_area.lr.x - written_subset.lr.x;

  if (size_above > 0.0) {
    write_stack->push_back(Area(old_area.ul.x,
                               old_area.ul.y,
                               old_area.lr.x,
                               written_subset.ul.y - 1));
  }

  if (size_below > 0.0) {
    write_stack->push_back(Area(old_area.ul.x,
                               written_subset.lr.y + 1,
                               old_area.lr.x,
                               old_area.lr.y));
  }

  if (size_left > 0.0) {
    write_stack->push_back(Area(old_area.ul.x,
                               written_subset.ul.y,
                               written_subset.lr.x - 1,
                               written_subset.lr.y));
  }

  if (size_right > 0.0) {
    write_stack->push_back(Area(written_subset.lr.x + 1,
                               written_subset.ul.y,
                               old_area.lr.x,
                               written_subset.lr.y));
  }

  return SP_None;
}

Area calculate_tile_intersection(PTIFF *tiff_file,
                                 Area subset) {
  const double tile_x_beginning = (subset.ul.x / tiff_file->block_x_size) * tiff_file->block_x_size;
  const double tile_y_beginning = (subset.ul.y / tiff_file->block_y_size) * tiff_file->block_y_size;
  const double tile_x_end = tile_x_beginning + tiff_file->block_x_size - 1;
  const double tile_y_end = tile_y_beginning + tiff_file->block_y_size - 1;
  const double subset_ul_x = std::max(tile_x_beginning, subset.ul.x);
  const double subset_ul_y = std::max(tile_y_beginning, subset.ul.y);
  const double subset_lr_x = std::min(tile_x_end, subset.lr.x);
  const double subset_lr_y = std::min(tile_y_end, subset.lr.y);

  return Area(subset_ul_x,
              subset_ul_y,
              subset_lr_x,
              subset_lr_y);  
}

SPTW_ERROR write_subset(PTIFF *tiff_file,
                        RasterChunk *chunk,
                        Area raster_subset) {
  const double tile_x_beginning = (raster_subset.ul.x / tiff_file->block_x_size) * tiff_file->block_x_size;
  const double tile_y_beginning = (raster_subset.ul.y / tiff_file->block_y_size) * tiff_file->block_y_size;
  const double tile_x_end = tile_x_beginning + tiff_file->block_x_size - 1;
  const double tile_y_end = tile_y_beginning + tiff_file->block_y_size - 1;

  MPI_Status status;
  int pixel_offset;
  MPI_Offset file_offset;

  int count = 0;
  Area ChunkArea;

  if (raster_subset.ul.x == tile_x_beginning 
      && raster_subset.lr.x == tile_x_end) {
    // Write whole subset in single call
    ChunkArea.ul = chunk->RasterToChunk(raster_subset.ul);
    ChunkArea.lr = chunk->RasterToChunk(raster_subset.lr);
    pixel_offset = (ChunkArea.ul.x + (ChunkArea.ul.y * chunk->column_count_))
        * tiff_file->band_type_size * chunk->band_count_;
    count = ((ChunkArea.lr.x - ChunkArea.ul.x + 1) * (ChunkArea.lr.y - ChunkArea.ul.y + 1))
        * tiff_file->band_type_size * chunk->band_count_;
    file_offset = chunk_to_file_offset(tiff_file,
                                       chunk,
                                       ChunkArea.ul.x + (ChunkArea.lr.y * chunk->column_count_));
    char *buffer = static_cast<char*>(malloc(count));
    for (int y = 0; y <= raster_subset.lr.y-raster_subset.ul.y; ++y) {
      int byte_row_size = chunk->column_count_*chunk->band_count_*tiff_file->band_type_size;
      int sub_row_size = (ChunkArea.lr.x - ChunkArea.ul.x) * tiff_file->band_type_size * chunk->band_count_;
      memcpy(buffer+(y*sub_row_size),
             static_cast<char*>(chunk->pixels_)+pixel_offset+(y*byte_row_size),
             sub_row_size);
    }
    MPI_File_write_at(tiff_file->fh,
                      file_offset,
                      buffer,
                      count,
                      MPI_BYTE,
                      &status);
    free(buffer);
  } else {
    // Write subset row by row
    const int row_count = raster_subset.lr.y - raster_subset.ul.y;
    for (int i = raster_subset.ul.y; i <= raster_subset.lr.y; ++i) {
      ChunkArea.ul = chunk->RasterToChunk(librasterblaster::Coordinate(raster_subset.ul.x, i, librasterblaster::UNDEF));
      ChunkArea.lr = chunk->RasterToChunk(librasterblaster::Coordinate(raster_subset.lr.x, i, librasterblaster::UNDEF));
      pixel_offset = (ChunkArea.ul.x + (ChunkArea.ul.y * chunk->column_count_))
          * tiff_file->band_type_size * chunk->band_count_;
      count = ((ChunkArea.lr.x - ChunkArea.ul.x + 1) * (ChunkArea.lr.y - ChunkArea.ul.y + 1))
          * tiff_file->band_type_size * chunk->band_count_;
      file_offset = chunk_to_file_offset(tiff_file,
                                         chunk,
                                         ChunkArea.ul.x + ChunkArea.ul.y * chunk->column_count_);

      MPI_File_write_at(tiff_file->fh,
                        file_offset,
                        (static_cast<char*>(chunk->pixels_)) + pixel_offset,
                        count,
                        MPI_BYTE,
                        &status);    
    }
  }
  
  return SP_None;
}

SPTW_ERROR write_rasterchunk(PTIFF *ptiff,
                             RasterChunk *chunk) {
  std::vector<Area> write_stack;
  Area write_area;
  write_area.ul = chunk->ChunkToRaster(librasterblaster::Coordinate(0.0, 0.0, librasterblaster::UNDEF));
  write_area.lr = chunk->ChunkToRaster(librasterblaster::Coordinate(chunk->column_count_-1,
                                                                    chunk->row_count_-1,
                                                                    librasterblaster::UNDEF));
  write_stack.push_back(write_area);

  while(!write_stack.empty()) {
    // Pop area needing write from top of stack
    Area top = write_stack.back();
    write_stack.pop_back();

    // Calculate subset of write area that is within the UL tile
    Area subset = calculate_tile_intersection(ptiff, top);

    // Fill the stack with any leftover areas
    fill_stack(&write_stack, top, subset);

    // Finally write the tile-bound subset
    write_subset(ptiff, chunk, subset);
  }
  SPTW_ERROR err = SP_None;
  unsigned char *pixels = NULL;
  return SP_None;
}
}
