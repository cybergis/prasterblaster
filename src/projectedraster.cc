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
// The ProjectedRaster class represents a raster with a location and a
// projection.
//


#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdlib>


#include "src/gctp_cpp/transformer.h"
#include "src/gctp_cpp/mercator.h"
#include "src/gctp_cpp/constants.h"
#include "src/gctp_cpp/utm.h"
#include "src/reprojection_tools.h"
#include "src/rasterchunk.h"
#include "src/sharedptr.h"
#include "src/projectedraster.h"

using std::string;

namespace librasterblaster {
ProjectedRaster::ProjectedRaster(string _filename) {
  OGRSpatialReference sr;

  filename_ = _filename;
  dataset_ = 0;
  row_count_ = column_count_ = -1;

  GDALAllRegister();

  ready_ = load_raster(filename_);

  return;
}

bool ProjectedRaster::CreateRaster(string _filename,
                                   int num_rows, int num_cols,
                                   GDALDataType pixel_type, double _pixel_size,
                                   int _band_count,
                                   shared_ptr<Projection> proj,
                                   double ulx, double uly) {
  bool status = make_raster(_filename,
                            num_cols,
                            num_rows,
                            _band_count,
                            ulx,
                            uly,
                            pixel_type,
                            proj,
                            _pixel_size);
  return status;
}

bool ProjectedRaster::CreateRaster(string _filename,
                                   shared_ptr<ProjectedRaster> input,
                                   shared_ptr<Projection> output_proj,
                                   double _pixel_size) {
  bool status = false;
  double ulx, uly;
  int num_cols, num_rows;

  Coordinate corner(input->ul_x(),
                    input->ul_y(),
                    UNDEF);

  Area out_area =  ProjectedMinbox(corner,
                                   input->srs(),
                                   input->pixel_size(),
                                   input->row_count(),
                                   input->column_count(),
                                   output_proj->wkt());

  ulx = out_area.ul.x;
  uly = out_area.ul.y;

  num_rows = static_cast<int>(ceil(out_area.ul.y - out_area.lr.y)
                              / _pixel_size);
  num_cols = static_cast<int>(ceil(out_area.lr.x - out_area.ul.x)
                              / _pixel_size);

  status = make_raster(_filename,
                      num_cols,
                      num_rows,
                      input->band_count_,
                      ulx,
                      uly,
                       input->pixel_type(),
                      output_proj,
                      _pixel_size);
  return status;
}

ProjectedRaster::~ProjectedRaster() {
  if (0 != dataset_) {
    GDALClose(static_cast<GDALDatasetH>(dataset_));
    dataset_ = 0;
  }
  return;
}

bool ProjectedRaster::ready() {
  if (projection_.get() != 0 &&
      (projection_->error() == 0) && ready_) {
    return true;
  }  else {
    return false;
  }
}

shared_ptr<Projection> ProjectedRaster::projection() {
  return projection_;
}


int ProjectedRaster::row_count() {
  return row_count_;
}

int ProjectedRaster::column_count() {
  return column_count_;
}

void ProjectedRaster::clamp_geo_coordinate(Coordinate *c) {
  if (c == 0) {
    return;
  }

  if (c->x < -180.0) {
    c->x = -179.999999;
  }

  if (c->x > 180.0) {
    c->x = 179.9999999;
  }

  if (c->y < -90.0) {
    c->y = -89.9999999;
  }

  if (c->y > 90.0) {
    c->y = 89.99999999;
  }

  return;
}

GDALDataType ProjectedRaster::pixel_type() {
  return pixel_type_;
}

int ProjectedRaster::bits_per_pixel()  {
  switch (pixel_type()) {
    case GDT_Byte:
      return 8;
    case GDT_Int16:
    case GDT_UInt16:
      return 16;
    case GDT_Int32:
    case GDT_UInt32:
    case GDT_Float32:
      return 32;
    case GDT_Float64:
      return 64;
    default:
      return -1;
  }
  return -1;
}

int ProjectedRaster::band_count() {
  if (dataset_ == 0) {
    return 1;
  } else {
    return dataset_->GetRasterCount();
  }

  return -1;
}

double ProjectedRaster::pixel_size() {
  return pixel_size_;
}

double ProjectedRaster::ul_x() {
  return ul_x_;
}

double ProjectedRaster::ul_y() {
  return ul_y_;
}

int ProjectedRaster::zone_number() {
  /*
    shared_ptr<UTM> utm_projection(projection);
	
    if (projection != 0) {
    if (projection->number() == _UTM) {
    return (utm_projection->zone());
    }
		
    }
  */
  return -1;
}

ProjDatum ProjectedRaster::datum() {
  if (projection_ != 0) {
    return projection_->datum();
  }
  return (ProjDatum)-1;
}

double* ProjectedRaster::gctp_parameters() {
  return gctp_parameters_;
}

bool ProjectedRaster::read_raster(int firstRow, int numRows, void *data) {
  bool success = false;

  if ((firstRow + numRows) > row_count()) {
    printf("Brrrrrangn\n\n\n");
    return false;
  }

  if (ready() && dataset_ != 0) {  // GTiff
    if (dataset_->RasterIO(GF_Read, 0, firstRow,
                          column_count_, numRows,
                          data, column_count_, numRows,
                          pixel_type_, band_count(),
                          NULL, 0, 0, 0) == CE_None) {
      success = true;
    }
  }
  return success;
}

bool ProjectedRaster::write_raster(int firstRow, int numRows, void *data) {
  bool success = false;

  if (firstRow < 0 || numRows +firstRow > row_count_ || data == 0) {
    fprintf(stderr, "Write Boink #1 %d %d\n", numRows + firstRow, row_count_);
    fflush(stderr);

    return false;
  }

  if (firstRow < 0 || firstRow >= dataset_->GetRasterYSize()) {
    fprintf(stderr, "Write boink #2\n");
    fflush(stderr);
    return false;
  }

  if ((firstRow + numRows) > dataset_->GetRasterYSize()) {
    fprintf(stderr, "Write boink #3\n");
    fflush(stderr);

    printf("Size: %d, Actual Size: %d\n", firstRow + numRows,
           dataset_->GetRasterYSize());
    return false;
  }
  // TODO(dmattli) Verify parameters!
  if (ready() && dataset_ != 0) {  // GTiff
    if (dataset_->RasterIO(GF_Write, 0, firstRow,
                          column_count_, numRows,
                          data, column_count_, numRows,
                          pixel_type_, band_count(),
                          NULL, 0, 0, 0) == CE_None) {
      success = true;
      dataset_->FlushCache();
    }
  }
  return success;
}

RasterChunk* ProjectedRaster::create_raster_chunk(Area area) {
RasterChunk *temp = create_allocated_raster_chunk(area);


  if (temp == NULL) {
    return NULL;
  }

  // Read area of raster
  if (ready() && dataset_ != 0) {
    if (dataset_->RasterIO(GF_Read,
                          area.ul.x,
                          area.ul.y,
                          temp->column_count_,
                          temp->row_count_,
                          temp->pixels_,
                          temp->column_count_,
                          temp->row_count_,
                          temp->pixel_type_,
                          band_count(),
                          NULL,
                          0, 0, 0) != CE_None) {
      // Cleanup and return error
      delete temp;
      return NULL;
    }
  } else {
    delete temp;
    return NULL;
  }

  return temp;
}

RasterChunk* ProjectedRaster::create_allocated_raster_chunk(Area area) {
  RasterChunk *temp = create_empty_raster_chunk(area);
  size_t buffer_size = (temp->row_count_ * temp->column_count_);
  temp->pixels_ = static_cast<unsigned char*>
      (calloc(buffer_size, GDALGetDataTypeSize(pixel_type())/8));

  if (temp->pixels_ == NULL) {
    fprintf(stderr, "Allocation error!\n");
    delete temp;
    return NULL;
  }

  return temp;
}

RasterChunk* ProjectedRaster::create_empty_raster_chunk(Area area) {
  RasterChunk *temp = new RasterChunk;

  temp->projection_ = shared_ptr<Projection>(projection());
  temp->raster_location_ = area.ul;
  temp->ul_projected_corner_ = Coordinate(ul_x_+(area.ul.x*pixel_size()),
                                          ul_y_-(area.ul.y*pixel_size()),
                                          UNDEF);
  temp->pixel_size_ = pixel_size();
  temp->row_count_ = area.lr.y - area.ul.y + 1;
  temp->column_count_ = area.lr.x - area.ul.x + 1;
  temp->pixel_type_ = pixel_type();
  temp->band_count_ = band_count_;
  temp->pixels_ = NULL;

  return temp;
}

bool ProjectedRaster::write_raster_chunk(RasterChunk *chunk) {
  // TODO(dmattli) Add some more checks

  if (dataset_->RasterIO(GF_Write,
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
    return false;
  }

  dataset_->FlushCache();

  return true;
}

bool ProjectedRaster::load_raster(string filename) {
  dataset_ = 0;
  row_count_ = column_count_ = -1;
  char *ref, **ugh;
  long projsys, zone, datum;
  double params[18] = {0.0};
  double *p;
  OGRSpatialReference sr;

  GDALAllRegister();

  dataset_ = static_cast<GDALDataset*>(GDALOpen(filename.c_str(), GA_Update));

  if (dataset_ == 0) {
    fprintf(stderr, "Error opening GDAL dataset\n");
    return false;
  }

  column_count_ = dataset_->GetRasterXSize();
  row_count_ = dataset_->GetRasterYSize();

  ref = strdup(dataset_->GetProjectionRef());
  if (ref == 0) {
    fprintf(stderr, "Error reading projection ref\n");
    return false;
  }

  double geo[6];
  dataset_->GetGeoTransform(geo);
  pixel_size_ = geo[1];
  ul_x_ = geo[0];
  ul_y_ = geo[3];
  pixel_type_ = (dataset_->GetRasterBand(1))->GetRasterDataType();
  band_count_ = dataset_->GetRasterCount();
  column_count_ = dataset_->GetRasterXSize();
  row_count_ = dataset_->GetRasterYSize();

  // Setup projection
  ugh = &ref;
  sr.importFromWkt(ugh);

  p = &(params[0]);
  sr.exportToUSGS(&projsys, &zone, &p, &datum);
  projection_ = shared_ptr<Projection>(
      Transformer::convertProjection(static_cast<ProjCode>(projsys)));

  if (projection_ == 0) {
    fprintf(stderr, "Error building projection, num %ld...\n", projsys);
    return false;
  }
  /*	shared_ptr<UTM> utm_projection(projection);
	if ((ProjCode)projsys == _UTM)
	(utm_ection->setZone(zone); */
  projection_->setParams(p);
  projection_->setDatum((ProjDatum)datum);

  projection_->setUnits(METER);
  ready_ = true;

  return true;
}

bool ProjectedRaster::make_raster(string _filename,
                                  int _cols,
                                  int _rows,
                                  int _band_count,
                                  double _ul_x,
                                  double _ul_y,
                                  GDALDataType _type,
                                  shared_ptr<Projection> _projection,
                                  double _pixel_size) {
  const char *format = "GTiff";
  GDALDriver *driver;
  char **options = 0;
  double geotransform[6] = { 444720, 30, 0, 3751320, 0, -30 };

  GDALAllRegister();

  driver = GetGDALDriverManager()->GetDriverByName(format);

  if (driver == NULL) {
    return false;
  }

  // Set options
  options = CSLSetNameValue(options, "INTERLEAVE", "PIXEL");
  options = CSLSetNameValue(options, "BIGTIFF", "YES");
  options = CSLSetNameValue(options, "TILED", "NO");
  options = CSLSetNameValue(options, "COMPRESS", "NONE");

  GDALDataset *_dataset = driver->Create(_filename.c_str(), _cols, _rows,
                                         _band_count, _type,
                                         options);

  if (_dataset == 0 || _projection == 0) {
    return false;
  }

  // Setup georeferencing
  OGRSpatialReference srs;

  geotransform[0] = _ul_x;
  geotransform[3] = _ul_y;
  geotransform[4] = geotransform[2] = 0.0;
  geotransform[1] = _pixel_size;
  geotransform[5] = -_pixel_size;
  _dataset->SetGeoTransform(geotransform);

  srs.importFromUSGS(_projection->number(),
                     0,
                     _projection->params(),
                     _projection->datum());

  char *wkt = NULL;
  srs.exportToWkt(&wkt);
  if (wkt != NULL) {
    _dataset->SetProjection(wkt);
    OGRFree(wkt);
  } else {
    fprintf(stderr,
            "\nERROR: Unsupported projection: %s\n",
            _projection->name().c_str());
  }

  GDALClose(_dataset);

  if (options != 0)
    CSLDestroy(options);

  return true;
}

string ProjectedRaster::srs() {
  OGRSpatialReference srs;
  string output_srs = "";
  char *temp;

  srs.importFromUSGS(projection_->number(),
                     0, projection_->params(),
                     projection_->datum());
  srs.exportToProj4(&temp);

  output_srs = temp;
  OGRFree(temp);

  return output_srs;
}
}
