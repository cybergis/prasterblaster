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
// The RasterChunk class represents an in-memory, georeferenced section of
// raster.
//
//

#include <gdal.h>

#include "src/rasterchunk.h"
#include "src/utils.h"

namespace librasterblaster {
RasterChunk* RasterChunk::CreateRasterChunk(GDALDataset *ds, Area chunk_area) {
  RasterChunk *temp = new RasterChunk;
  double gt[6];
  ds->GetGeoTransform(gt);
  ds->GetGeoTransform(temp->geotransform_);

  temp->projection_ = ds->GetProjectionRef();
  temp->raster_location_ = chunk_area.ul;
  temp->ul_projected_corner_ = Coordinate(gt[0]+(chunk_area.ul.x*gt[1]),
                                          gt[3]-(chunk_area.ul.y*gt[1]),
                                          UNDEF);
  temp->pixel_size_ = gt[1];
  temp->row_count_ = chunk_area.lr.y - chunk_area.ul.y + 1;
  temp->column_count_ = chunk_area.lr.x - chunk_area.ul.x + 1;
  temp->pixel_type_ = ds->GetRasterBand(1)->GetRasterDataType();
  temp->band_count_ = ds->GetRasterCount();
  temp->pixels_ = NULL;

  size_t buffer_size = (temp->row_count_ * temp->column_count_);
  temp->pixels_ = static_cast<unsigned char*>
      (calloc(buffer_size, GDALGetDataTypeSize(temp->pixel_type_)/8));

  if (temp->pixels_ == NULL) {
    fprintf(stderr, "Allocation error!\n");
    delete temp;
    return NULL;
  }

  return temp;
}

PRB_ERROR RasterChunk::ReadRasterChunk(GDALDataset *ds, RasterChunk *chunk) {
  // Read area of raster

  if (ds == NULL) {
    return PRB_BADARG;
  }

  if (ds->RasterIO(GF_Read,
                   chunk->raster_location_.x,
                   chunk->raster_location_.y,
                   chunk->column_count_,
                   chunk->row_count_,
                   chunk->pixels_,
                   chunk->column_count_,
                   chunk->row_count_,
                   chunk->pixel_type_,
                   chunk->band_count_,
                   NULL,
                   0, 0, 0) != CE_None) {
    return PRB_IOERROR;
  }

  return PRB_NOERROR;
}

PRB_ERROR RasterChunk::WriteRasterChunk(GDALDataset *ds, RasterChunk *chunk) {
  if (ds->RasterIO(GF_Write,
                   chunk->raster_location_.x,
                   chunk->raster_location_.y,
                   chunk->column_count_,
                   chunk->row_count_,
                   chunk->pixels_,
                   chunk->column_count_,
                   chunk->row_count_,
                   chunk->pixel_type_,
                   chunk->band_count_,
                   NULL, 0, 0, 0) != CE_None) {
    // Error!
    fprintf(stderr, "Error writing RasterChunk %p\n", chunk->pixels_);
    return PRB_IOERROR;
  }

  ds->FlushCache();

  return PRB_NOERROR;
}

}
