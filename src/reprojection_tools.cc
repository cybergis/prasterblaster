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
// Helper functions to create and manipulate projections and projected rasters.
//
//

#include "src/reprojection_tools.h"

#include <float.h>
#include <stdint.h>

#include <ogr_api.h>
#include <ogr_spatialref.h>
#include <gdal.h>
#include <gdal_priv.h>

#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/coordinate.h"

#include "src/rastercoordtransformer.h"
#include "src/projectedraster.h"
#include "src/quadtree.h"
#include "src/resampler.h"
#include "src/sharedptr.h"
#include "src/utils.h"


namespace librasterblaster {
PRB_ERROR CreateOutputRaster(GDALDataset *in,
                             string output_filename,
                             double output_pixel_size,
                             string output_srs) {
  OGRSpatialReference in_srs;
  OGRSpatialReference out_srs;
  OGRErr err;

  err = in_srs.SetFromUserInput(in->GetProjectionRef());
  if (err != OGRERR_NONE) {
    return PRB_BADARG;
  }

  out_srs.SetFromUserInput(output_srs.c_str());
  if (err != OGRERR_NONE) {
    return PRB_BADARG;
  }

  // Determine output raster size by calculating the projected coordinate minbox
  double in_transform[6];
  in->GetGeoTransform(in_transform);
  Coordinate ul(in_transform[0], in_transform[3], UNDEF);
  char *srs_str = NULL;
  in_srs.exportToProj4(&srs_str);
  Area out_area = ProjectedMinbox(ul,
                                  srs_str,
                                  in_transform[1],
                                  in->GetRasterYSize(),
                                  in->GetRasterXSize(),
                                  output_srs);
  CPLFree(srs_str);

  int64_t num_rows = static_cast<int64_t>(ceil(out_area.ul.y - out_area.lr.y)
                              / in_transform[1]);
  int64_t num_cols = static_cast<int64_t>(ceil(out_area.lr.x - out_area.ul.x)
                              / in_transform[1]);

  GDALAllRegister();
  GDALDriver *driver = GetGDALDriverManager()->GetDriverByName("GTiff");

  if (driver == NULL) {
    fprintf(stderr, "Error opening GTiff driver.\n");
    return PRB_BADARG;
  }

  // Set driver options
  char **options = NULL;
  options = CSLSetNameValue(options, "INTERLEAVE", "PIXEL");
  options = CSLSetNameValue(options, "BIGTIFF", "YES");
  options = CSLSetNameValue(options, "TILED", "NO");
  options = CSLSetNameValue(options, "COMPRESS", "NONE");


  GDALDataset *output =
      driver->Create(output_filename.c_str(),
                     num_cols,
                     num_rows,
                     in->GetRasterCount(),
                     in->GetRasterBand(1)->GetRasterDataType(),
                     options);

  if (output == NULL) {
    fprintf(stderr, "driver->Create call failed.\n");
    return PRB_BADARG;
  }

  // Copy ColorTable
  if (in->GetRasterCount() > 0) {
    GDALColorTable *ct = in->GetRasterBand(1)->GetColorTable();
    output->GetRasterBand(1)->SetColorTable(ct);
  }

  // Setup georeferencing
  double out_t[6] = { out_area.ul.x,
                      in_transform[1],
                      0.0,
                      out_area.ul.y,
                      0.0,
                      in_transform[5] };

  output->SetGeoTransform(out_t);
  OGRSpatialReference out_sr;
  char *wkt = NULL;
  out_sr.SetFromUserInput(output_srs.c_str());
  out_sr.exportToWkt(&wkt);
  output->SetProjection(wkt);

  OGRFree(wkt);
  CSLDestroy(options);
  GDALClose(output);

  return PRB_NOERROR;
}

PRB_ERROR CreateSampleOutput(shared_ptr<ProjectedRaster> input,
                             string output_filename,
                             string output_srs,
                             int output_size) {
  shared_ptr<Projection> in_proj = input->projection();
  shared_ptr<Projection> out_proj;
  OGRSpatialReference srs;

  OGRErr err = srs.importFromProj4(output_srs.c_str());

  if (err != OGRERR_NONE) {
    fprintf(stderr, "Error parsing projection!\n");
    return PRB_PROJERROR;
  }

  int64_t proj_code, datum_code, zone;
  double *params = NULL;

  srs.exportToUSGS(&proj_code, &zone, &params, &datum_code);

  out_proj = shared_ptr<Projection>(
      Transformer::convertProjection(static_cast<ProjCode>(proj_code)));

  if (!out_proj) {
    return PRB_PROJERROR;
  }

  out_proj->setUnits(in_proj->units());
  out_proj->setDatum(in_proj->datum());
  out_proj->setParams(params);

  OGRFree(params);

  Coordinate ul(input->ul_x(), input->ul_y(), UNDEF);
  Area parea = ProjectedMinbox(ul,
                               in_proj->wkt(),
                               input->pixel_size(),
                               input->row_count(),
                               input->column_count(),
                               output_srs);


  int xsize = (parea.lr.x - parea.ul.x) / output_size;
  int ysize = (parea.ul.y - parea.lr.y) / output_size;
  double pixel_size = (parea.lr.x - parea.ul.x) / xsize;

  if (ysize > xsize) {
    xsize = ysize;
    pixel_size = (parea.ul.y - parea.lr.y) / xsize;
  }

  bool result = ProjectedRaster::CreateRaster(output_filename,
                                              input,
                                              out_proj,
                                              pixel_size);


  if (result) {
    return PRB_NOERROR;
  } else {
    return PRB_IOERROR;
  }
}

Projection* ProjectionFactory(string output_srs) {
  Projection *out_proj;
  OGRSpatialReference srs;
  char *wkt = NULL;
  char *tmp = NULL;

  OGRErr err = srs.importFromProj4(output_srs.c_str());

  if (err != OGRERR_NONE) {
    wkt = strdup(output_srs.c_str());
    tmp = wkt;
    err = srs.importFromWkt(&tmp);
    free(wkt);
    if (err != OGRERR_NONE) {
      fprintf(stderr, "Error parsing projection!\n");
      return NULL;
    }
  }


  int64_t proj_code, datum_code, zone;
  double *params = NULL;

  srs.exportToUSGS(&proj_code, &zone, &params, &datum_code);

  out_proj = Transformer::convertProjection(static_cast<ProjCode>(proj_code));

  if (!out_proj) {
    return NULL;
  }

  out_proj->setDatum(static_cast<ProjDatum>(datum_code));
  out_proj->setParams(params);
  out_proj->setUnits(METER);

  OGRFree(params);

  return out_proj;
}

// For use by PartitionBySize
int simplerandom(int i) {
  return std::rand()%i;
}

std::vector<Area> PartitionBySize(int rank,
                                  int process_count,
                                  int row_count,
                                  int column_count,
                                  int maximum_partition_size) {
  // Seed the psuedorandom generator. This _must_ be the same on all processes.
  std::srand(42);

  QuadTree qt(row_count,
              column_count,
              maximum_partition_size);

  vector<Area> leaves = qt.collectLeaves();
  vector<Area> partitions;
  size_t partition_count = leaves.size();
  size_t partitions_per_proc = partition_count / process_count;

  // Now we shuffle the partitions for load balancing.
  // All processes should generate the same shuffle.
  std::random_shuffle(leaves.begin(), leaves.end(), simplerandom);

  if (partitions_per_proc == 0) {
    // process_count < partition_count
    if (static_cast<size_t>(rank) < partition_count) {
      partitions.push_back(leaves.at(rank));
    }
  } else {
    // process_count >= partition_count
    size_t first_index = rank * partitions_per_proc;
    size_t last_index = (rank+1) * partitions_per_proc - 1;

    if (rank == process_count - 1) {
      last_index = leaves.size() - 1;
    }

    for (size_t i = first_index; i <= last_index; ++i) {
      partitions.push_back(leaves.at(i));
    }
  }

  return partitions;
}


void SearchAndUpdate(Area input_area,
                     shared_ptr<Projection> input_projection,
                     shared_ptr<Projection> output_projection,
                     double input_ulx,
                     double input_uly,
                     double input_pixel_size,
                     Area *output_area) {
  Coordinate input_coord;
  Coordinate temp;
  OGRSpatialReference input_sr, output_sr;
  char *input_wkt = strdup(input_projection->wkt().c_str());
  char *output_wkt = strdup(output_projection->wkt().c_str());
  char *intemp = input_wkt;
  char *outtemp = output_wkt;
  input_sr.importFromWkt(&intemp);
  output_sr.importFromWkt(&outtemp);
  OGRCoordinateTransformation *t =
      OGRCreateCoordinateTransformation(&input_sr,
                                        &output_sr);

  if (t == NULL) {
          output_area->ul.x = -1.0;
          output_area->ul.y = -1.0;
          return;
  }

  for (int64_t x = input_area.ul.x; x <= input_area.lr.x; ++x) {
    for (int64_t y = input_area.ul.y; y >= input_area.lr.y; --y) {
      input_coord.x = x * input_pixel_size + input_ulx;
      input_coord.y = input_uly - (y * input_pixel_size);

      /*
      input_projection->inverse(input_coord.x,
                                input_coord.y, &temp.x, &temp.y);
      output_projection->forward(temp.x, temp.y, &temp.x, &temp.y);
      */
      t->Transform(1, &input_coord.x, &input_coord.y);
      temp = input_coord;

      if (temp.x  < output_area->ul.x) {
        output_area->ul.x = temp.x;
      }
      if (temp.y > output_area->ul.y) {
        output_area->ul.y = temp.y;
      }
      if (temp.x > output_area->lr.x)
        output_area->lr.x = temp.x;
      if (temp.y < output_area->lr.y)
        output_area->lr.y = temp.y;
    }
  }

  if (input_wkt != NULL) {
    free(input_wkt);
  }


  if (output_wkt != NULL) {
    free(output_wkt);
  }

  OCTDestroyCoordinateTransformation(t);
  return;
}

Area ProjectedMinbox(Coordinate input_ul_corner,
                     string input_srs,
                     double input_pixel_size,
                     int input_row_count,
                     int input_column_count,
                     string output_srs) {
  // Input area, projected coordinates
  Area ia;
  // Projected Area
  Area output_area;
  shared_ptr<Projection> input_proj(ProjectionFactory(input_srs));
  shared_ptr<Projection> output_projection(ProjectionFactory(output_srs));
  const int buffer = 2;

  if (input_proj.get() == NULL || output_projection.get() == NULL) {
    return Area(-1, -1, -1, -1);
  }

  output_area.ul.x = output_area.lr.y = DBL_MAX;
  output_area.ul.y = output_area.lr.x = -DBL_MAX;

  // Check the top of the raster
  ia.ul.x = 0;
  ia.lr.x = input_column_count - 1;
  ia.ul.y = buffer;
  ia.lr.y = 0;

  SearchAndUpdate(ia,
                  input_proj,
                  output_projection,
                  input_ul_corner.x,
                  input_ul_corner.y,
                  input_pixel_size,
                  &output_area);

  // Check the bottom of the raster
  ia.ul.x = 0;
  ia.lr.x = input_column_count - 1;
  ia.ul.y = input_row_count - 1;
  ia.lr.y = input_row_count - 1 - buffer;

  SearchAndUpdate(ia,
                  input_proj,
                  output_projection,
                  input_ul_corner.x,
                  input_ul_corner.y,
                  input_pixel_size,
                  &output_area);

  // Check Left
  ia.ul.x = 0;
  ia.lr.x = buffer;
  ia.ul.y = input_row_count - 1;
  ia.lr.y = 0;

  SearchAndUpdate(ia,
                  input_proj,
                  output_projection,
                  input_ul_corner.x,
                  input_ul_corner.y,
                  input_pixel_size,
                  &output_area);

  // Check right
  ia.ul.x = input_column_count - buffer - 1;
  ia.lr.x = input_column_count - 1;
  ia.ul.y = input_row_count - 1;
  ia.lr.y = 0;

  SearchAndUpdate(ia,
                  input_proj,
                  output_projection,
                  input_ul_corner.x,
                  input_ul_corner.y,
                  input_pixel_size,
                  &output_area);

  return output_area;
}
Area RasterMinbox(GDALDataset *source,
                  GDALDataset *destination,
                  Area destination_raster_area) {
  double s_gt[6];
  double d_gt[6];
  source->GetGeoTransform(s_gt);
  destination->GetGeoTransform(d_gt);

  Coordinate s_ul(s_gt[0], s_gt[3], UNDEF);
  Coordinate d_ul(d_gt[0], d_gt[3], UNDEF);

  string s_srs(source->GetProjectionRef());
  string d_srs(destination->GetProjectionRef());

  return RasterMinbox2(s_srs,
                       s_ul,
                       s_gt[1],
                       source->GetRasterYSize(),
                       source->GetRasterXSize(),
                       d_srs,
                       d_ul,
                       d_gt[1],
                       destination->GetRasterYSize(),
                       destination->GetRasterXSize(),
                       destination_raster_area);
}

Area RasterMinbox2(string source_projection,
                  Coordinate source_ul,
                  double source_pixel_size,
                  int source_row_count,
                  int source_column_count,
                  string destination_projection,
                  Coordinate destination_ul,
                  double destination_pixel_size,
                  int destination_row_count,
                  int destination_column_count,
                  Area destination_raster_area) {
  Area source_area;
  Coordinate c;
  RasterCoordTransformer rt(source_projection,
                            source_ul,
                            source_pixel_size,
                            source_row_count,
                            source_column_count,
                            destination_projection,
                            destination_ul,
                            destination_pixel_size);
  if (rt.ready() == false) {
    return Area(-1.0, -1.0, -1.0, -1.0);
  }

  Area temp;
  int buffer = 5;

  source_area.ul.x = source_area.ul.y = DBL_MAX;
  source_area.lr.y = source_area.lr.x = -DBL_MAX;

  int step = 1;

  for (int x = destination_raster_area.ul.x;
       x <= destination_raster_area.lr.x; x+=step) {
    for (int y = destination_raster_area.ul.y;
         y <= destination_raster_area.lr.y; ++y) {
      c.x = x;
      c.y = y;

      temp = rt.Transform(c);

      if (temp.ul.x == -1) {
        continue;
      }

      // Check that calculated minbox in within destination raster space.
      if ((temp.ul.x < -0.01) || (temp.ul.x > destination_column_count)
          || (temp.ul.y < 0.0) || (temp.ul.y > destination_row_count)
          || (temp.lr.x > destination_column_count - 1) || (temp.lr.x < 0.0)
          || (temp.lr.y > destination_row_count -1) || (temp.lr.y < 0.0)) {
        temp.ul.x = -1.0;
        continue;
        printf("\n\nSearch Area: %f %f %f %f\n",
               destination_raster_area.ul.x,
               destination_raster_area.ul.y,
               destination_raster_area.lr.x,
               destination_raster_area.lr.y);
        printf("Source UL: %f %f\n", source_ul.x, source_ul.y);
        printf("Destin UL: %f %f\n", destination_ul.x, destination_ul.y);
        printf("Source raster size, columns: %d, rows %d\n",
               destination_column_count, destination_row_count);
        printf("Source: %d %d\n", x, y);
        printf("Outside rasterspace: %f %f %f %f\n",
               temp.ul.x, temp.ul.y, temp.lr.x, temp.lr.y);
      }

      if (temp.lr.x > source_area.lr.x) {
        source_area.lr.x = temp.lr.x;
      }

      if (temp.ul.x > source_area.lr.x) {
        source_area.lr.x = temp.ul.x;
      }

      if (temp.ul.x < source_area.ul.x) {
        source_area.ul.x = temp.ul.x;
      }

      if (temp.lr.x < source_area.ul.x) {
        source_area.ul.x = temp.lr.x;
      }

      if (temp.ul.y < source_area.ul.y) {
        source_area.ul.y = temp.ul.y;
      }

      if (temp.lr.y > source_area.lr.y) {
        source_area.lr.y = temp.lr.y;
      }
    }
  }

  // Check whether entire area is out of the projected space.
  if ((source_area.ul.x == DBL_MAX) || (source_area.ul.y == DBL_MAX)
      || (source_area.lr.x == -DBL_MAX) || (source_area.lr.y == -DBL_MAX)) {
    source_area.ul.x = -1.0;
    source_area.lr.x = -1.0;
    source_area.ul.y = -1.0;
    source_area.lr.y = -1.0;
    return source_area;
  }


  source_area.ul.x = floor(source_area.ul.x);
  source_area.ul.y = floor(source_area.ul.y);
  source_area.lr.x = ceil(source_area.lr.x);
  source_area.lr.y = ceil(source_area.lr.y);

  if (source_area.lr.x > destination_column_count - 1) {
    source_area.lr.x = destination_column_count - 1;
  }

  if (source_area.lr.y > destination_row_count - 1) {
    source_area.lr.y = destination_row_count - 1;
  }

  if (source_area.lr.y < source_area.ul.y) {
    source_area.lr.y = source_area.ul.y;
  }

  if (source_area.lr.x < source_area.ul.x) {
    source_area.lr.x = source_area.ul.x;
  }

  return source_area;
}

/**
 * \brief This function takes two RasterChunk pointers and performs
 *        reprojection and resampling
 * \param source Pointer to the RasterChunk to reproject from
 * \param destination Pointer to the RasterChunk to reproject to
 * \param fillvalue std::string with the fillvalue
 * \param resampler The resampler that should be used
 *
 * @return Returns a bool indicating success or failure.
 */
bool ReprojectChunk(RasterChunk *source,
                    RasterChunk *destination,
                    string fillvalue,
                    RESAMPLER resampler) {
  if (source->pixel_type_ != destination->pixel_type_) {
    fprintf(stderr, "Source and destination chunks have different types!\n");
    return false;
  }

  double fvalue = strtod(fillvalue.c_str(), NULL);

  switch (source->pixel_type_) {
    case GDT_Byte:
      GEN_RESAMPLER_CASES(uint8_t);
    case GDT_UInt16:
      GEN_RESAMPLER_CASES(uint16_t);
    case GDT_Int16:
      GEN_RESAMPLER_CASES(int16_t);
    case GDT_UInt32:
      GEN_RESAMPLER_CASES(uint32_t);
    case GDT_Int32:
      GEN_RESAMPLER_CASES(int32_t);
    case GDT_Float32:
      GEN_RESAMPLER_CASES(float);
    case GDT_Float64:
      GEN_RESAMPLER_CASES(double);
    case GDT_CInt16:
    case GDT_CInt32:
    case GDT_CFloat32:
    case GDT_CFloat64:
    default:
      fprintf(stderr, "Invalid type in ReprojectChunk!\n");
      return false;
      break;
  }
  return true;
}
}
