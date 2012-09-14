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

#include <algorithm>
#include <float.h>
#include <vector>
#include <stdint.h>

#include <ogr_api.h>
#include <ogr_spatialref.h>
#include <gdal.h>
#include <gdal_priv.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/coordinate.h"

#include "rastercoordtransformer.h"o
#include "projectedraster.h"
#include "quadtree.h"
#include "resampler.h"
#include "sharedptr.h"
#include "utils.h"

namespace librasterblaster {
PRB_ERROR CreateOutputRaster(shared_ptr<ProjectedRaster> in,
                             string output_filename,
                             double output_pixel_size,
                             string output_srs)
{
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
		
  out_proj = shared_ptr<Projection>(Transformer::convertProjection(static_cast<ProjCode>(proj_code)));
	
  if (!out_proj) {
    return PRB_PROJERROR;
  }

  out_proj->setUnits(in_proj->units());
  out_proj->setDatum(in_proj->datum());
  out_proj->setParams(params);

  OGRFree(params);

  bool result  = ProjectedRaster::CreateRaster(output_filename,
                                               in,
                                               shared_ptr<Projection>(out_proj->copy()),
                                               in->pixel_type() ,
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
                             int output_size)
{
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
		
  out_proj = shared_ptr<Projection>(Transformer::convertProjection(static_cast<ProjCode>(proj_code)));
	
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
                                              input->pixel_type(),
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
	
  OGRErr err = srs.importFromProj4(output_srs.c_str());
	
  if (err != OGRERR_NONE) {
    fprintf(stderr, "Error parsing projection!\n");
    return NULL;
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

  OGRFree(params);

  return out_proj;

}

std::vector<Area> PartitionByCount(shared_ptr<ProjectedRaster> source,
				   int partition_count) {
  Area noval(-1, -1, -1, -1);
  int64_t area_size = static_cast<int64_t>(source->column_count()) * static_cast<int64_t>(source->row_count());
  int64_t partition_area = area_size / partition_count;
	
  if (partition_area == 0) {
    partition_area = 1;
  }

  QuadTree tree(source->row_count(), source->column_count(), partition_area);

  return tree.collectLeaves();
}


void SearchAndUpdate(Area input_area,
                     shared_ptr<Projection> input_projection,
                     shared_ptr<Projection> output_projection,
                     double input_ulx,
                     double input_uly,
                     double input_pixel_size,
                     Area *output_area)
{
  Coordinate input_coord;
  Coordinate temp;

  for (int64_t x = input_area.ul.x; x <= input_area.lr.x; ++x) {
    for (int64_t y = input_area.ul.y; y >= input_area.lr.y; --y) {
      input_coord.x = x * input_pixel_size + input_ulx;
      input_coord.y = input_uly - (y * input_pixel_size);
			
      input_projection->inverse(input_coord.x, 
                                input_coord.y, &temp.x, &temp.y);
      output_projection->forward(temp.x, temp.y, &temp.x, &temp.y);
			
      if (temp.x  < output_area->ul.x) 
        output_area->ul.x = temp.x;
      if (temp.y > output_area->ul.y)
        output_area->ul.y = temp.y;
      if (temp.x > output_area->lr.x) 
        output_area->lr.x = temp.x;
      if (temp.y < output_area->lr.y)
        output_area->lr.y = temp.y;
    }
  }

  return;
}

Area ProjectedMinbox(Coordinate input_ul_corner,
                     string input_srs,
                     double input_pixel_size,
                     int input_row_count,
                     int input_column_count,
                     string output_srs) {
  Area ia; // Input area, projected coorsdinates
  Area output_area; // Projected Area
  shared_ptr<Projection> input_proj(ProjectionFactory(input_srs));
  shared_ptr<Projection> output_projection(ProjectionFactory(output_srs));
  const int buffer = 2;
  
  if (input_proj.get() == NULL || output_projection.get()) {
    return Area(-1, -1, -1, -1);
  }

  output_area.ul.x = output_area.lr.y = DBL_MAX;
  output_area.ul.y = output_area.lr.x = -DBL_MAX;
  
  // Check the top of the raster
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

  // Check the bottom of the raster
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
  ia.ul.x = input_column_count - 1;
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
		  Area destination_raster_area)
{
	
  Area source_area;
  Coordinate c(destination->ul_x(), destination->ul_y(), UNDEF);
  shared_ptr<Projection> dproj = destination->projection();
  RasterCoordTransformer rt(source, dproj,
                            c, destination->pixel_size());
  Area temp;


  source_area.ul.x = source_area.ul.y = DBL_MAX;
  source_area.lr.y = source_area.lr.x = -DBL_MAX;
	
  for (int x = destination_raster_area.ul.x; x <= destination_raster_area.lr.x; ++x) {
    for (int y = destination_raster_area.ul.y; y <= destination_raster_area.lr.y; ++y) {

      c.x = x;
      c.y = y;

      temp = rt.Transform(c);
		
      if (temp.ul.x == -1) {
        continue;
      }

      if (temp.lr.x > source_area.lr.x) {
        source_area.lr.x = temp.lr.x;
      }
			
      if (temp.ul.x < source_area.ul.x) {
        source_area.ul.x = temp.ul.x;
      }

      if (temp.ul.y < source_area.ul.y) {
        source_area.ul.y = temp.ul.y;
      }

      if (temp.lr.y > source_area.lr.y) {
        source_area.lr.y = temp.lr.y;
      }
    }
  }

  if (source_area.ul.x < 0)
    source_area.ul.x = 0;
  if (source_area.ul.y < 0)
    source_area.ul.y = 0;
	
  source_area.ul.x = floor(source_area.ul.x);
  source_area.ul.y = floor(source_area.ul.y);
  source_area.lr.x = ceil(source_area.lr.x);
  source_area.lr.y = ceil(source_area.lr.y);

  if (source_area.lr.x > destination->column_count() - 1) {
    source_area.lr.x = destination->column_count() - 1;
  }

  if (source_area.lr.y > destination->row_count() - 1) {
    source_area.lr.y = destination->row_count() - 1;
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
                                               in->pixel_type(),
                                               output_pixel_size);
  if (result) {
    return PRB_NOERROR;
  } else {
    return PRB_NOERROR;
  }
}

bool ReprojectChunk(RasterChunk *source, 
                    RasterChunk *destination, 
                    string fillvalue, 
                    string resampler_name) {
  if (source->pixel_type_ != destination->pixel_type_) {
    fprintf(stderr, "Source and destination chunks have different types!\n");
    return false;
  }

  double fvalue = strtod(fillvalue.c_str(), NULL);

  std::transform(resampler_name.begin(), resampler_name.end(), resampler_name.begin(), ::tolower);

  RESAMPLER resampler = MIN;

  if (resampler_name == "min") {
    resampler = MIN;
  } else if (resampler_name == "max") {
    resampler = MAX;
  } else if (resampler_name == "nearest") {
    resampler = NEAREST;
  }
	
  switch (source->pixel_type_) 
  {
    case GDT_Byte:
      switch (resampler){
        case MIN:
          return ReprojectChunkType<unsigned char>(source, destination, static_cast<uint8_t>(fvalue), &(Min<unsigned char>));
          break;
        case MAX:
          return ReprojectChunkType<unsigned char>(source, destination, static_cast<uint8_t>(fvalue), &(Max<unsigned char>));
        case NEAREST:
        default:
          return ReprojectChunkType<unsigned char>(source, destination, static_cast<uint8_t>(fvalue), NULL);
			
      }
      break;
    case GDT_UInt16:
      return ReprojectChunkType<unsigned short>(source, destination, static_cast<uint16_t>(fvalue), &(Max<unsigned short>));
      break;
    default:
      fprintf(stderr, "Invalid type in ReprojectChunk!\n");
      return false;
      break;

  }
  return true;
}
}
