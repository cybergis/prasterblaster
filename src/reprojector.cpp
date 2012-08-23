/*!
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
 * 
 *
 */

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <stdint.h>

#include <ogr_api.h>
#include <ogr_spatialref.h>
#include <gdal.h>
#include <gdal_priv.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/coordinate.h"

#include "reprojector.hh"
#include "resampler.hh"
#include "sharedptr.hh"
#include "quadtree.hh"


RasterCoordTransformer::RasterCoordTransformer(shared_ptr<ProjectedRaster> source,
					       shared_ptr<ProjectedRaster> _dest)
					       
{
	src_proj = shared_ptr<Projection>(source->getProjection());
	source_pixel_size_ = source->getPixelSize();
	source_ul_ = Coordinate(source->ul_x, source->ul_y, UNDEF);
	dest_proj = shared_ptr<Projection>(_dest->getProjection());
	destination_pixel_size_ = _dest->getPixelSize();
	destination_ul_ = Coordinate(_dest->ul_x, _dest->ul_y, UNDEF);
	
	return;
}

RasterCoordTransformer::RasterCoordTransformer(shared_ptr<ProjectedRaster> source,
					       shared_ptr<Projection> destination_projection,
					       Coordinate destination_ul,
					       double destination_pixel_size) 

{
	destination_pixel_size_ = destination_pixel_size;
	destination_ul_ = destination_ul;
	source_ul_ = Coordinate(source->ul_x, source->ul_y, UNDEF);
	source_pixel_size_ = source->getPixelSize();

	dest_proj = destination_projection;
	src_proj =  shared_ptr<Projection>(source->getProjection());

	return;
}
					       
   
RasterCoordTransformer::RasterCoordTransformer(shared_ptr<Projection> source_projection,
					       Coordinate source_ul,
					       double source_pixel_size,
					       shared_ptr<Projection> destination_projection,
					       Coordinate destination_ul,
					       double destination_pixel_size)
{
	src_proj = source_projection;
	source_ul_ = source_ul;
	source_pixel_size_ = source_pixel_size;
	dest_proj = destination_projection;
	destination_ul_ = destination_ul;
	destination_pixel_size_ = destination_pixel_size;

}

RasterCoordTransformer::~RasterCoordTransformer()
{
	return;
}

Area RasterCoordTransformer::Transform(Coordinate source)
{

	Area value;
	Coordinate temp1, temp2;
	
	temp1.x = temp1.y = 0.0;

	value.ul = temp1;
	value.lr = temp1;

	temp1.x = (static_cast<double>(source.x) * source_pixel_size_) + source_ul_.x;
	temp1.y = (static_cast<double>(source.y) * source_pixel_size_) - source_ul_.y;

	src_proj->inverse(temp1.x, temp1.y, &temp2.x, &temp2.y);
	src_proj->forward(temp2.x, temp2.y, &temp2.x, &temp2.y);

	if (fabs(temp1.x - temp2.x) > 0.01) {
		// Point is outside defined projection area, return no-value
		value.ul.x = -1.0;
		value.lr.x = -1.0;
		return value;
	}

	temp1.x = (static_cast<double>(source.x) * source_pixel_size_) + source_ul_.x;
	temp1.y = (static_cast<double>(source.y) * source_pixel_size_) - source_ul_.y;
	temp2 = temp1;

	// Now we are going to assign temp1 as the UL of our pixel and
	// temp2 as LR
//	temp1.x -= src->pixel_size/2;
//	temp1.y += src->pixel_size/2;
	temp2.x += source_pixel_size_/2;
	temp2.y -= source_pixel_size_/2;

	src_proj->inverse(temp1.x, temp1.y, &temp1.x, &temp1.y);
	dest_proj-> forward(temp1.x, temp1.y, &temp1.x, &temp1.y);
	src_proj->inverse(temp2.x, temp2.y, &temp2.x, &temp2.y);
	dest_proj-> forward(temp2.x, temp2.y, &temp2.x, &temp2.y);
	// temp1/temp2 now contain coords to input projection
	// Now convert to points in the raster coordinate space.
	temp1.x -= destination_ul_.x;
	temp1.y += destination_ul_.y;
	temp1.x /= destination_pixel_size_;
	temp1.y /= destination_pixel_size_;
	temp2.x -= destination_ul_.x;
	temp2.y += destination_ul_.y;
	temp2.x /= destination_pixel_size_;
	temp2.y /= destination_pixel_size_;

	value.ul = temp1;
	value.lr = temp2;
	
	return value;
}

PRB_ERROR CreateOutputRaster(shared_ptr<ProjectedRaster> in,
			string output_filename,
			double output_pixel_size,
			string output_srs)
{
	shared_ptr<Projection> in_proj = shared_ptr<Projection>(in->getProjection());
	shared_ptr<Projection> out_proj;
	OGRSpatialReference srs; 
	
	OGRErr err = srs.importFromProj4(output_srs.c_str());
	
	if (err != OGRERR_NONE) {
		fprintf(stderr, "Error parsing projection!\n");
		return PROJ_ERROR;
	}

	long proj_code, datum_code, zone;
	double *params = NULL;
		
	srs.exportToUSGS(&proj_code, &zone, &params, &datum_code);
		
	out_proj = shared_ptr<Projection>(Transformer::convertProjection(static_cast<ProjCode>(proj_code)));
	
	if (!out_proj) {
		return PROJ_ERROR;
	}

	out_proj->setUnits(in_proj->units());
	out_proj->setDatum(in_proj->datum());
	out_proj->setParams(params);

	OGRFree(params);

	bool result  = ProjectedRaster::CreateRaster(output_filename,
						     in,
						     shared_ptr<Projection>(out_proj->copy()),
						     in->type ,
						     output_pixel_size);
	if (result) {
		return NO_ERROR;
	} else {
		return FILE_ERROR;
	}


}

PRB_ERROR CreateSampleOutput(shared_ptr<ProjectedRaster> input,
			string output_filename,
			string output_srs, 
			int output_size)
{
	shared_ptr<Projection> in_proj = input->getProjection();
	shared_ptr<Projection> out_proj;
	OGRSpatialReference srs; 
	
	OGRErr err = srs.importFromProj4(output_srs.c_str());
	
	if (err != OGRERR_NONE) {
		fprintf(stderr, "Error parsing projection!\n");
		return PROJ_ERROR;
	}

	long proj_code, datum_code, zone;
	double *params = NULL;
		
	srs.exportToUSGS(&proj_code, &zone, &params, &datum_code);
		
	out_proj = shared_ptr<Projection>(Transformer::convertProjection(static_cast<ProjCode>(proj_code)));
	
	if (!out_proj) {
		return PROJ_ERROR;
	}

	out_proj->setUnits(in_proj->units());
	out_proj->setDatum(in_proj->datum());
	out_proj->setParams(params);

	OGRFree(params);
	
	Area parea = ProjectedMinbox(input,
				     out_proj);


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
						    input->type,
						    pixel_size);


	if (result) {
		return NO_ERROR;
	} else {
		return FILE_ERROR;
	}
}

PRB_ERROR WriteRasterChunk(string output_filename,
			   RasterChunk::RasterChunk *chunk)
{
	Coordinate t;
	t.x = chunk->raster_location_.x;
	t.y = chunk->raster_location_.y;
	/*
	ProjectedRaster::CreateRaster(output_filename,
				      chunk->row_count_,
				      chunk->column_count_,
				      chunk->pixel_type_,
				      chunk->pixel_size_,
				      chunk->band_count_);

	*/
	chunk->raster_location_.x = 0;
	chunk->raster_location_.y = 0;
	ProjectedRaster *out = new ProjectedRaster(output_filename);
	out->writeRasterChunk(chunk);
	chunk->raster_location_.x = t.x;
	chunk->raster_location_.y = t.y;
	

	

	return NO_ERROR;
}

std::vector<Area> PartitionByCount(shared_ptr<ProjectedRaster> source,
				   int partition_count)
{
	Area noval(-1, -1, -1, -1);
	int64_t area_size = static_cast<int64_t>(source->getColCount()) * static_cast<int64_t>(source->getRowCount());
	int64_t partition_area = area_size / partition_count;
	
	if (partition_area == 0) {
		partition_area = 1;
	}

	QuadTree tree(source->getRowCount(), source->getColCount(), partition_area);

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
		     
		     

Area ProjectedMinbox(shared_ptr<ProjectedRaster> input,
		     shared_ptr<Projection> output_projection)
{
	Area ia; // Input area, projected coordinates
	Area output_area; // Projected Area
	shared_ptr<Projection> input_proj(input->getProjection());
	const int buffer = 2;
	
	output_area.ul.x = output_area.lr.y = DBL_MAX;
	output_area.ul.y = output_area.lr.x = -DBL_MAX;
	
	// Check the top of the raster
	ia.ul.x = 0;
	ia.lr.x = input->getColCount() - 1;
	ia.ul.y = input->getRowCount() - 1;
	ia.lr.y = input->getRowCount() - 1 - buffer;
	
	SearchAndUpdate(ia,
			input_proj,
			output_projection,
			input->ul_x,
			input->ul_y,
			input->getPixelSize(),
			&output_area);

	// Check the bottom of the raster
	ia.ul.x = 0;
	ia.lr.x = input->getColCount() - 1;
	ia.ul.y = buffer;
	ia.lr.y = 0;

	SearchAndUpdate(ia,
			input_proj,
			output_projection,
			input->ul_x,
			input->ul_y,
			input->getPixelSize(),
			&output_area);

	// Check Left
	ia.ul.x = 0;
	ia.lr.x = buffer;
	ia.ul.y = input->getRowCount() - 1;
	ia.lr.y = 0;

	SearchAndUpdate(ia,
			input_proj,
			output_projection,
			input->ul_x,
			input->ul_y,
			input->getPixelSize(),
			&output_area);



	// Check right
	ia.ul.x = input->getColCount() - 1;
	ia.lr.x = input->getColCount() - 1;
	ia.ul.y = input->getRowCount() - 1;
	ia.lr.y = 0;
	
	SearchAndUpdate(ia,
			input_proj,
			output_projection,
			input->ul_x,
			input->ul_y,
			input->getPixelSize(),
			&output_area);

	return output_area;
}

Area RasterMinbox(shared_ptr<ProjectedRaster> source,
		  shared_ptr<ProjectedRaster> destination,
		  Area destination_raster_area)
{
	
	Area source_area;
	RasterCoordTransformer rt(source, destination->getProjection(), 
				  Coordinate(destination->ul_x, destination->ul_y, UNDEF), destination->getPixelSize());
	Area temp;
	Coordinate c;

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

	if (source_area.lr.x > destination->getColCount() - 1) {
		source_area.lr.x = destination->getColCount() - 1;
	}

	if (source_area.lr.y > destination->getRowCount() - 1) {
		source_area.lr.y = destination->getRowCount() - 1;
	}

	return source_area;
}

bool ReprojectChunk(RasterChunk::RasterChunk *source, RasterChunk::RasterChunk *destination, string fillvalue, string resampler_name)
{
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
			return ReprojectChunkType<unsigned char>(source, destination, static_cast<uint8_t>(fvalue), &(Resampler::Min<unsigned char>));
			break;
		case MAX:
			return ReprojectChunkType<unsigned char>(source, destination, static_cast<uint8_t>(fvalue), &(Resampler::Max<unsigned char>));
		case NEAREST:
		default:
			return ReprojectChunkType<unsigned char>(source, destination, static_cast<uint8_t>(fvalue), NULL);
			
		}
		break;
	case GDT_UInt16:
		return ReprojectChunkType<unsigned short>(source, destination, static_cast<uint16_t>(fvalue), &(Resampler::Max<unsigned short>));
		break;
	default:
		fprintf(stderr, "Invalid type in ReprojectChunk!\n");
		return false;
		break;

	}
	return true;
}

void updateMinbox(double x, double y, Area *minbox) 
{
	if (x > -DBL_MAX && x < minbox->ul.x) 
		minbox->ul.x = x;
	if (x < DBL_MAX && x > minbox->lr.x) 
		minbox->lr.x = x;
	if (y > -DBL_MAX && y < minbox->lr.y) 
		minbox->lr.y = y;
	if (y < DBL_MAX && y > minbox->ul.y) 
		minbox->ul.y = y;
	
	return;
}

Area FindRasterArea(shared_ptr<ProjectedRaster> source_raster,
		    shared_ptr<ProjectedRaster> dest_raster,
		    int first_row,
		    int last_row)
{
	RasterCoordTransformer trans(source_raster, dest_raster);

	Area dest_area;
	Coordinate temp;
	Area temp_area;
	
	dest_area.ul.x = dest_area.lr.y = -DBL_MAX;
	dest_area.lr.x = dest_area.ul.y = DBL_MAX;
	
	for (int row = first_row; row <= last_row; ++row) {
		for (int col = 0; col < source_raster->getColCount(); ++col) {
			temp.x = col;
			temp.y = row;

			temp_area = trans.Transform(temp);
			if (temp_area.ul.x == -1.0) {
				continue;
			}
			
			updateMinbox(temp_area.ul.x, temp_area.ul.y, &dest_area);
			updateMinbox(temp_area.lr.x, temp_area.lr.y, &dest_area);
		}
	}

	return dest_area;
}

void clampGeoCoordinate(Coordinate *c)
{
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
		    

