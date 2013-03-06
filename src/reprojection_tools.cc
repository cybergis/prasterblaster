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
PRB_ERROR CreateOutputRaster(shared_ptr<ProjectedRaster> in,
                             string output_filename,
                             double output_pixel_size,
                             string output_srs) {
  shared_ptr<Projection> in_proj = shared_ptr<Projection>(in->projection());
  shared_ptr<Projection> out_proj;
  OGRSpatialReference srs;

  OGRErr err = srs.importFromProj4(output_srs.c_str());

  if (err != OGRERR_NONE) {
    fprintf(stderr, "Error parsing projection: %s\n", output_srs.c_str());
    return PRB_PROJERROR;
  }

  long proj_code, datum_code, zone;
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

  bool result  = ProjectedRaster::CreateRaster(output_filename,
                                               in,
                                               out_proj,
                                               output_pixel_size);
  if (result) {
    return PRB_NOERROR;
  } else {
    return PRB_IOERROR;
  }
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

  long proj_code, datum_code, zone;
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


  long proj_code, datum_code, zone;
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

std::vector<Area> RowPartition(int rank,
                               int process_count,
                               int row_count,
                               int column_count,
                               long partition_size) {
  int rows_per_partition = partition_size / column_count;
  int partitions_per_row = column_count / partition_size;
  int leftover_per_row = column_count % partition_size;
  int partition_count = 1;

  if (rows_per_partition > 0) {
    partition_count = row_count / rows_per_partition;
  }
  vector<Area> partitions;

  if (partitions_per_row < 1) {
    partitions_per_row = 0;
    partition_size = partition_size / column_count;
    leftover_per_row = 0;

    if (rows_per_partition > row_count) {
      rows_per_partition = 1;
    }
  }

  if (rows_per_partition > 0) {
    for (int i=0; i<partition_count-1; ++i) {
      partitions.push_back(Area(0,
                                i*rows_per_partition,
                                column_count-1,
                                (i+1)*rows_per_partition-1));
    }
    // Final partition, may have > rows_per_partition
    partitions.push_back(Area(0,
                              partition_count * rows_per_partition,
                              column_count-1,
                              row_count-1));
    return partitions;
    
  }
  
  vector<int> rowpart_sizes(partitions_per_row, partition_size);
  vector<int> first_column(partitions_per_row, 0);

  for (int i = 0; i < leftover_per_row; ++i) {
    // Set the partition sizes
    int j = i % rowpart_sizes.size();
    rowpart_sizes.at(j) = rowpart_sizes.at(j) + 1;
  }

  // Set the column offsets
  for (int i = 1; i < first_column.size(); ++i) {
    first_column.at(i) = first_column.at(i-1) + rowpart_sizes.at(i-1);
  }

  partition_count = rowpart_sizes.size() * row_count;

  for (int r = rank; r < partition_count; r += process_count) {
    int row = r / rowpart_sizes.size();
    int partcol = r % rowpart_sizes.size();

    partitions.push_back(Area(first_column.at(partcol),
                              row,
                              first_column.at(partcol)
                              + rowpart_sizes.at(partcol) - 1,
                              row));
  }

  std::reverse(partitions.begin(), partitions.end());
  return partitions;
}

std::vector<Area> PartitionBySize(int rank,
                                  int process_count,
                                  int row_count,
                                  int column_count,
                                  int partition_size,
                                  int maximum_height,
                                  int maximum_width) {
  if (maximum_height == -1) {
    maximum_height = row_count + 1;
  }

  if (maximum_width == -1) {
    maximum_width = column_count + 1;
  }

  QuadTree qt(row_count,
              column_count,
              partition_size);

  vector<Area> leaves = qt.collectLeaves();
  vector<Area> partitions;
  size_t partition_count = leaves.size();
  size_t partitions_per_proc = partition_count / process_count;

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
  OGRCoordinateTransformation *t = OGRCreateCoordinateTransformation(&input_sr,
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
        //printf("New UL_Y: %f %f, %f %f\n", temp.x, temp.y, x * input_pixel_size + input_ulx, input_uly - (y * input_pixel_size));

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


Area RasterMinbox(shared_ptr<ProjectedRaster> source,
                  shared_ptr<ProjectedRaster> destination,
                  Area destination_raster_area) {
  Coordinate s_ul(source->ul_x(),
                  source->ul_y(),
                  UNDEF);

  Coordinate d_ul(destination->ul_x(),
                  destination->ul_y(),
                  UNDEF);
  shared_ptr<Projection> s_proj = source->projection();
  shared_ptr<Projection> d_proj = destination->projection();
  return RasterMinbox(s_proj,
                      s_ul,
                      source->pixel_size(),
                      source->row_count(),
                      source->column_count(),
                      d_proj,
                      d_ul,
                      destination->pixel_size(),
                      destination->row_count(),
                      destination->column_count(),
                      destination_raster_area);
}
Area RasterMinbox(shared_ptr<Projection> source_projection,
                  Coordinate source_ul,
                  double source_pixel_size,
                  int source_row_count,
                  int source_column_count,
                  shared_ptr<Projection> destination_projection,
                  Coordinate destination_ul,
                  double destination_pixel_size,
                  int destination_row_count,
                  int destination_column_count,
                  Area destination_raster_area) {
  Area source_area;
  Coordinate c;
  shared_ptr<Projection> dproj = destination_projection;
  RasterCoordTransformer rt(source_projection,
                            source_ul,
                            source_pixel_size,
                            source_row_count,
                            source_column_count,
                            destination_projection,
                            destination_ul,
                            destination_pixel_size);
  Area temp;
  int buffer = 5;

  if (dproj->errorOccured() == true) {
    fprintf(stderr, "Error with destination projection in RasterMinbox\n");
    source_area.ul.x = -1.0;
    source_area.lr.x = -1.0;
    return source_area;
  }

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
//      printf("\t Source point: %d %d\n", x, y);
//      printf("\t\t Calculated point: %f %f %f %f\n", temp.ul.x, temp.ul.y, temp.lr.x, temp.lr.y);
      if (temp.ul.x == -1) {
        continue;
      }
      
      // Check that calculated minbox in within destination raster space.
      if ((temp.ul.x < -0.01) || (temp.ul.x > destination_column_count)
          || (temp.ul.y < 0.0) || (temp.ul.y > destination_row_count)
          || (temp.lr.x > destination_column_count - 1) || (temp.lr.x < 0.0)
          || (temp.lr.y > destination_row_count -1) || (temp.lr.y < 0.0)) {
        printf("Source raster size, rows: %f, columns %f\n", destination_raster_area.lr.x, destination_raster_area.lr.y);
        printf("Source: %d %d\n", x, y);
        printf("Outside rasterspace: %f %f %f %f\n", temp.ul.x, temp.ul.y, temp.lr.x, temp.lr.y);

        
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
    double t = source_area.lr.y;
    source_area.lr.y = source_area.ul.y;
  }

  if (source_area.lr.x < source_area.ul.x) {
    source_area.lr.x = source_area.ul.x;
  }

  return source_area;
}

PRB_ERROR CreateOutputRaster2(shared_ptr<ProjectedRaster> in,
                             string output_filename,
                             double output_pixel_size,
                             string output_srs) {
  shared_ptr<Projection> in_proj = shared_ptr<Projection>(in->projection());
  shared_ptr<Projection> out_proj;
  OGRSpatialReference srs;

  OGRErr err = srs.importFromProj4(output_srs.c_str());

  if (err != OGRERR_NONE) {
    fprintf(stderr, "Error parsing projection!\n");
    return PRB_PROJERROR;
  }

  long proj_code, datum_code, zone;
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

  bool result  = ProjectedRaster::CreateRaster(output_filename,
                                               in,
                                               shared_ptr<Projection>(
                                                   out_proj->copy()),
                                               output_pixel_size);
  if (result) {
    return PRB_NOERROR;
  } else {
    return PRB_NOERROR;
  }
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
      switch (resampler) {
        case MIN:
          return ReprojectChunkType<unsigned char>(source,
                                                   destination,
                                                   static_cast<uint8_t>(fvalue),
                                                   &(Min<unsigned char>));
          break;
        case MAX:
          return ReprojectChunkType<unsigned char>(source,
                                                   destination,
                                                   static_cast<uint8_t>(fvalue),
                                                   &(Max<uint8_t>));
    case NEAREST:
    default:
          return ReprojectChunkType<unsigned char>(source,
                                                   destination,
                                                   static_cast<uint8_t>(fvalue),
                                                   NULL);
      }
      break;
    case GDT_UInt16:
      return ReprojectChunkType<uint16_t>(source,
                                                destination,
                                                static_cast<uint16_t>(fvalue),
                                                &(Max<uint16_t>));
      break;
    default:
      fprintf(stderr, "Invalid type in ReprojectChunk!\n");
      return false;
      break;
  }
  return true;
}
}
