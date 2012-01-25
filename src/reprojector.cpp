/*!
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE
 * This work was produced as a part of the official duties of a
 * federal employee and is in the public domain.

 * @section DESCRIPTION
 *
 *
 */
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <cstdint>
#include <sstream>
#include <vector>

#include <mpi.h>
#include <gdal.h>
#include <gdal_priv.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/coordinate.h"

#include "rastercache.hh"
#include "reprojector.hh"
#include "resampler.hh"


using std::shared_ptr;


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

	temp1.x = ((double)source.x * source_pixel_size_) + source_ul_.x;
	temp1.y = ((double)source.y * source_pixel_size_) - source_ul_.y;
	
	src_proj->inverse(temp1.x, temp1.y, &temp2.x, &temp2.y);
	src_proj->forward(temp2.x, temp2.y, &temp2.x, &temp2.y);

	if (fabs(temp1.x - temp2.x) > 0.01) {
		value.ul.x = std::numeric_limits<double>::quiet_NaN();
		value.lr.x = value.ul.x;
	}

	temp1.x = ((double)source.x * source_pixel_size_) + source_ul_.x;
	temp1.y = ((double)source.y * source_pixel_size_) - source_ul_.y;
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

vector<Area> PartitionByCount(shared_ptr<ProjectedRaster> source,
			       int partition_count)
{
	vector<Area> partitions(partition_count);
	int chunk_size = (source->getColCount() * source->getRowCount()) / partition_count;
	Area temp;
	
	for (int i = 0; i < (int)partitions.size(); ++i) {
		temp.ul.x = (chunk_size * i) % source->getColCount();
		temp.ul.y = (chunk_size * i) / source->getColCount();
		temp.lr.x = temp.ul.x + chunk_size % source->getColCount();

		if (temp.lr.x == 0.0) {  // Wrap column values around 
			temp.lr.x = source->getColCount() - 1;
		} else {
			--temp.lr.x;
		}
		
		temp.lr.y = temp.ul.y + chunk_size / source->getColCount() - 1;

		partitions.at(i) = temp;
	}

	partitions.back().lr.x = source->getColCount() - 1;
	partitions.back().lr.y = source->getRowCount() - 1;

	return partitions;


}

Area FindOutputArea(shared_ptr<ProjectedRaster> input,
		      shared_ptr<Projection> output_projection,
		      double output_pixel_size)
{
	Area geoarea; // Projected Area
	shared_ptr<Projection> input_proj(input->getProjection());
	
	geoarea.ul.x = geoarea.lr.y = DBL_MAX;
	geoarea.ul.y = geoarea.lr.x = -DBL_MAX;
	
	for (int x = 0; x < input->getColCount(); ++x) {
		for (int y = 0; y < input->getRowCount(); ++y) {
			Coordinate input_coord;
			Coordinate temp;

			input_coord.x = x * input->getPixelSize() + input->ul_x;
			input_coord.y = y * input->getPixelSize() * input->ul_y;
			
			input_proj->inverse(input_coord.x, input_coord.y, &temp.x, &temp.y);
			output_projection->forward(temp.x, temp.y, &temp.x, &temp.y);

			if (temp.x  < geoarea.ul.x) 
				geoarea.ul.x = temp.x;
			if (temp.y > geoarea.ul.y)
				geoarea.ul.y = temp.y;
			if (temp.x > geoarea.lr.x) 
				geoarea.lr.x = temp.x;
			if (temp.y < geoarea.lr.y)
				geoarea.lr.y = temp.y;
					
		}

	}

	
	
	return geoarea;
}

Area MapDestinationAreatoSource(shared_ptr<ProjectedRaster> source,
				shared_ptr<Projection> destination_projection,
				Area destination_area,
				double destination_pixel_size)
{
	
	Area source_area;
	RasterCoordTransformer rt(source, destination_projection, destination_area.ul, destination_pixel_size);
	Area temp;
	Coordinate c;
        int first(0), last(0);

	source_area.lr.x = 0;
	source_area.lr.y = 0;
	
	source_area.ul.x = source->getColCount();
	source_area.ul.y = source->getRowCount();

	for (int x = destination_area.ul.x; x < destination_area.lr.x; ++x) {
		for (int y = destination_area.ul.y; y <= destination_area.lr.y; ++y) {

			c.x = x;
			c.y = y;

			temp = rt.Transform(c);
                
			if (c.x > source_area.lr.x) {
				source_area.lr.x = c.x;
			}
			
			if (c.x < source_area.ul.x) {
				source_area.ul.x = c.x;
			}

			if (c.y < source_area.ul.y) {
				source_area.ul.y = c.y;
			}

			if (c.y > source_area.lr.y) {
				source_area.lr.y = c.y;
			}
		}
	}

	return source_area;
}

bool ParallelReprojection(shared_ptr<ProjectedRaster> source, shared_ptr<ProjectedRaster> destination, 
			  int rank, int process_count)
{
	vector<Area> parts = PartitionByCount(destination, process_count);
	int parts_per_process = parts.size() / process_count;
	int leftover = parts.size() % process_count;

	

}

bool ReprojectChunk(RasterChunk::RasterChunk source, RasterChunk::RasterChunk destination)
{
	if (source.pixel_type_ != destination.pixel_type_) {
		return false;
	}


	switch (source.pixel_type_) 
	{
	case GDT_Byte:
		ReprojectChunkType<unsigned char>(source, destination);
		break;
	case GDT_UInt16:
		ReprojectChunkType<uint16_t>(source, destination);
		break;
	default:
		return false;
		break;

	}
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
			try {
				temp_area = trans.Transform(temp);
			} catch(std::string) {
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
		    

Area FindRasterExtent(shared_ptr<ProjectedRaster> raster,
		      Area geographical_area)
{

	Area rasterArea;
	double pixel_size = raster->pixel_size;
	shared_ptr<Projection> projection(raster->getProjection());
	Area projArea;// = FindProjectedExtent(projection, geographical_area, raster->pixel_size);
	int rows = (int)((projArea.lr.x - projArea.ul.x) / pixel_size);
	int cols = (int)((projArea.ul.y - projArea.lr.y) / pixel_size);
	int first_row = (int)((raster->ul_x - projArea.ul.x) / pixel_size);
	int first_col = 0;

	rasterArea.ul.y = first_row;
	rasterArea.ul.x = first_col;
	rasterArea.lr.y = first_row + rows;
	rasterArea.lr.x =  cols;

	return rasterArea;

}
