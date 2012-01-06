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
#include <stdexcept>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <set>
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
	src = source;
	dest = _dest;
	src_proj = shared_ptr<Projection>(source->getProjection());
	dest_proj = shared_ptr<Projection>(dest->getProjection());
	
	return;
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

	temp1.x = ((double)source.x * src->pixel_size) + src->ul_x;
	temp1.y = ((double)source.y * src->pixel_size) - src->ul_y;
	
	src_proj->inverse(temp1.x, temp1.y, &temp2.x, &temp2.y);
	src_proj->forward(temp2.x, temp2.y, &temp2.x, &temp2.y);

	if (fabs(temp1.x - temp2.x) > 0.01) {
		value.ul.x = std::numeric_limits<double>::quiet_NaN();
		value.lr.x = value.ul.x;
	}

	temp1.x = ((double)source.x * src->pixel_size) + src->ul_x;
	temp1.y = ((double)source.y * src->pixel_size) - src->ul_y;
	temp2 = temp1;

	// Now we are going to assign temp1 as the UL of our pixel and
	// temp2 as LR
//	temp1.x -= src->pixel_size/2;
//	temp1.y += src->pixel_size/2;
	temp2.x += src->pixel_size/2;
	temp2.y -= src->pixel_size/2;

	src_proj->inverse(temp1.x, temp1.y, &temp1.x, &temp1.y);
	dest_proj-> forward(temp1.x, temp1.y, &temp1.x, &temp1.y);
	src_proj->inverse(temp2.x, temp2.y, &temp2.x, &temp2.y);
	dest_proj-> forward(temp2.x, temp2.y, &temp2.x, &temp2.y);
	// temp1/temp2 now contain coords to input projection
	// Now convert to points in the raster coordinate space.
	temp1.x -= dest->ul_x;
	temp1.y += dest->ul_y;
	temp1.x /= dest->pixel_size;
	temp1.y /= dest->pixel_size;
	temp2.x -= dest->ul_x;
	temp2.y += dest->ul_y;
	temp2.x /= dest->pixel_size;
	temp2.y /= dest->pixel_size;

	value.ul = temp1;
	value.lr = temp2;
	
	return value;
}

vector<Interval> ParititionByCount(shared_ptr<ProjectedRaster> source,
				   int partition_count)
{
	vector<Interval> partitions(partition_count);
	int source_size = (source->getColCount() * source->getRowCount()) / partition_count;
	int last_size = (source->getColCount() * source->getRowCount()) % partition_count;
	Interval temp;
	
	if (last_size == 0) {
		last_size = source_size;
	}

	for (int i = 0; i < partitions.size()-1; ++i) {

		temp.first_index = i * partition_count;
		temp.last_index = ((i+1) * partition_count) - 1;
		partitions.at(i) = (temp);
	}
	// Set last interval
	temp.first_index = (partitions.size()-1) * partition_count;
	temp.last_index = (source->getColCount() * source->getRowCount()) - 1;
	partitions.back() = temp;

	return partitions;
}

Interval FindSourceInterval(shared_ptr<ProjectedRaster> source,
			    shared_ptr<ProjectedRaster> destination,
			    Interval dest_area)
{
	Interval source_interval(source->getRowCount()*source->getColCount(), 0);
	RasterCoordTransformer rt(source, destination);
	Area temp;
	Coordinate c;
        int first(0), last(0);

	for (int i = dest_area.first_index; i <= dest_area.last_index; ++i) {
		c.x = i % destination->getColCount();
		c.y = i / destination->getColCount();

		temp = rt.Transform(c);
                
                first = (int)temp.ul.x * temp.ul.y * source->getColCount();
                last = (int) temp.lr.x * temp.lr.y * source->getColCount();

                if (first < source_interval.first_index) {
                        source_interval.first_index = first;
                }

                if (last > source_interval.last_index) {
                        source_interval.last_index = last;
                }
	}

	return source_interval;
}


bool ReprojectChunk(RasterChunk::RasterChunk source, RasterChunk::RasterChunk destination)
{
        shared_ptr<Projection> outproj, inproj;
        Coordinate temp1, temp2;
	std::vector<char> inraster, outraster;
	set<long> overflow_rows;
/*	

        outproj = output->getProjection();
        inproj = input->getProjection();
	
        // Setup raster vectors
	outraster.resize(chunk.rowCount * output->getColCount() * output->bitsPerPixel()/8);
	printf("Process: %u of %d Output raster chunk %ld %ld\n",  
	       rank, numprocs, chunk.firstIndex, chunk.firstIndex+chunk.rowCount); 
	Area pixelArea;
	int count = 0;
	int total = 0;
	RasterCoordTransformer rt(output, input);        

	for (int chunk_y = 0; chunk_y < chunk.rowCount; ++chunk_y)  {
		for (int chunk_x = 0; chunk_x < output->getColCount(); ++chunk_x) {

			// Determine location of equivalent input pixel
			temp1.x = chunk_x; 
			temp1.y = chunk_y + chunk.firstIndex;
			try {
				pixelArea = rt.Transform(temp1);   /// Now makes sense.
			} catch (std::runtime_error) {
				continue;
			}
			temp1 = pixelArea.ul;
			temp2 = pixelArea.lr;
			if (temp1.x < 0)
				temp1.x = 0;
			// temp1&2 are now scaled to input raster coords, now resample!
			// But does the rectangle defined by temp1 and temp2 actually
			// contain any points? If not use nearest-neighbor...
			long center_index = (long)((temp1.x + temp2.x)/2);
			center_index += (long)((temp1.y + temp2.y)/2) * input->getColCount();
			long pixel_width = (long)fabs((temp2.x - temp1.x)/input->getPixelSize());
			long pixel_height = (long)(fabs(temp1.y - temp2.y)/input->getPixelSize());

			if (pixel_width == 0) {
 				pixel_width = 1;
			}
			if (pixel_height == 0) {
				pixel_height = 1;
			}
			std::vector<long> index_array(pixel_width*pixel_height, 0L);

			temp1.x = floor(temp1.x);
			temp1.y = floor(temp1.y);

			long index = 0;
			unsigned char val1 = 0;

			index = (unsigned int)temp1.x;
			index += ((unsigned int)temp1.y) * input->getColCount();
			try{
				val1 = sourceCache.at(index);
			} catch (out_of_range &e) {
				printf("%s\n", e.what());

			}
			try {
				outraster.at(chunk_x+(chunk_y*output->getColCount())) = val1;
				++total;
			} catch (out_of_range &e) {
				printf("Chunk_x %d, Chunk_y %d, up to %d %ld\n", chunk_x, chunk_y,
				       output->getColCount(), chunk.rowCount);
				outraster.at(chunk_x+(chunk_y*output->getColCount())) = 255;
				++count;

			}

		
		}
	}

	// Write output raster
	output->writeRaster(chunk.firstIndex, chunk.rowCount, &(outraster[0]));
*/	
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
	Area projArea = FindProjectedExtent(projection, geographical_area, raster->pixel_size);
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
