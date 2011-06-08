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
#include <set>
#include <sstream>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <mpi.h>
#include <gdal.h>
#include <gdal_priv.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/coordinate.h"

#include "rastercache.hh"
#include "reprojector.hh"
#include "resampler.hh"


using boost::shared_ptr;


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
		std::ostringstream s;
		s << "Point Overlaps! " << source.x << ", " << source.x;
		// Overlap detected!
		throw std::runtime_error(s.str());
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

Chunker::Chunker(shared_ptr<ProjectedRaster> source, 
		shared_ptr<ProjectedRaster> destination)
{
	source_raster = source;
	dest_raster = destination;
	return;
}

void Chunker::clampGeoCoordinate(Coordinate *c)
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


vector<ChunkExtent> Chunker::getChunksByCount(int process_count, int process_index) 
{
	std::vector<ChunkExtent> chunks;
	long chunk_size = dest_raster->getRowCount() / process_count;
	Area geominbox, projminbox;

	if (process_count <= 0) {
		throw std::invalid_argument("Invalid process count");
	}

	if (dest_raster->getRowCount() <= 0) {
		throw std::runtime_error("Destination raster has no rows");
	}

	for (int i=0; i < process_count-1; ++i) {
		chunks.push_back(ChunkExtent(i*chunk_size,
					     ((i+1)*chunk_size-1)));
	}

	// Last chunk includes any leftover
	long last_row = dest_raster->getRowCount() - 1;
	long last_chunk_size = last_row - (process_count - 1) + 1;

	chunks.push_back(ChunkExtent((process_count-1)*chunk_size, last_row));


	return chunks;
}
	

Reprojector::Reprojector(shared_ptr<ProjectedRaster> _input, shared_ptr<ProjectedRaster> _output,
			 int processCount, int rank) 
        :
        numprocs(processCount), rank(rank), input(_input), output(_output)
{
	numchunks = numprocs * 10;
        maxx = maxy = 0;
        minx = miny = 1e+37;
        resampler = &resampler::nearest_neighbor<unsigned char>;


        return;
};

Reprojector::~Reprojector()
{

        return;
}

void Reprojector::parallelReproject()
{
	Chunker chunker(input, output);
	printf("Finding chunks\n");
	vector<ChunkExtent> chunks = chunker.getChunksByCount(numprocs, rank);
	
	printf("Reprojecting chunks...\n");
	for (int i = 0; i < chunks.size(); ++i) {
		reprojectChunk(chunks[i]);
	}


}

void Reprojector::reprojectChunk(ChunkExtent chunk)
{
        shared_ptr<Projection> outproj, inproj;
        Coordinate temp1, temp2;
	std::vector<char> inraster, outraster;
	set<long> overflow_rows;
	RasterCache<unsigned char> sourceCache(input);
	

        outproj = output->getProjection();
        inproj = input->getProjection();
	
        // Setup raster vectors
	outraster.resize(chunk.rowCount * output->getColCount() * output->bitsPerPixel()/8);
 
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

			//! TODO: Check for overlap!
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
			try {
				index = (unsigned int)temp1.x;
				index += ((unsigned int)temp1.y) * input->getColCount();

				outraster.at(chunk_x+(chunk_y*output->getColCount())) = sourceCache.at(index);
				++total;
			} catch (out_of_range &e) {
				overflow_rows.insert(temp1.y);
				outraster.at(chunk_x+(chunk_y*output->getColCount())) = 255;
				++count;

			}

		
		}
	}
	printf("done!\n");
	printf("%d total, %d pixels off! %f percent\n", total, count, (double)count/(double)total*100);
	for (set<long>::iterator i=overflow_rows.begin(); i != overflow_rows.end(); ++i) {
//		printf("Tried to access %d\n", *i);
	}
	// Write output raster
	output->writeRaster(chunk.firstIndex, chunk.rowCount, &(outraster[0]));
	
	return;
}

void Reprojector::reproject()
{
	parallelReproject();

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
		    

Area FindGeographicalExtent(shared_ptr<Projection> projection, 
			    Coordinate ul_point,
			    int row_count,
			    int col_count,
			    double pixel_size)
{
	Area geoarea;
	double lon, lat;
	double x, y;
	int check_count = 5;
	Coordinate temp;
	
	geoarea.ul.x = DBL_MAX;
	geoarea.ul.y = -DBL_MAX;
	geoarea.lr.x = -DBL_MAX;
	geoarea.lr.y = DBL_MAX;
	lon = lat = x = y = 0.0;

	if (row_count < check_count) {
		check_count = row_count;
	}

	// Check top
	for (int row = 0; row < check_count; row += 5) {
		for (int col = 0; col < col_count; ++col) {
			x = ul_point.x + (pixel_size * col);
			y = ul_point.y - (pixel_size * row);
			projection->inverse(x, y, &lon, &lat);
			temp.x = lon;
			temp.y = lat;
			projection->forward(temp.x, temp.y, &temp.x, &temp.y);
			if (fabs(temp.x - x) > 0.001) {
				continue;
			} else {
				updateMinbox(lon, lat, &geoarea);
			}
		}
	}
	// Check bottom
	for (int row = row_count-1; row > row_count-1-check_count; row -= 5) {
		for (int col = 0; col < col_count; ++col) {
			x = ul_point.x + (pixel_size * col);
			y = ul_point.y - (pixel_size * row);
			projection->inverse(x, y, &lon, &lat);
			temp.x = lon;
			temp.y = lat;
			projection->forward(temp.x, temp.y, &temp.x, &temp.y);
			if (fabs(temp.x - x) > 0.0001) {
				continue;
			} else {
				updateMinbox(lon, lat, &geoarea);
			}
		}
	}

	clampGeoCoordinate(&geoarea.ul);
	clampGeoCoordinate(&geoarea.lr);

	return geoarea;

}

Area FindProjectedExtent(shared_ptr<Projection> projection,
			 Area geographical_area,
			 double pixel_size)
{
	Area projarea;
	Coordinate projected;

	projarea.ul.x = DBL_MAX;
	projarea.ul.y = -DBL_MAX;
	projarea.lr.x = -DBL_MAX;
	projarea.lr.y = DBL_MAX;

	// Check top of map
	for (double lon = geographical_area.ul.x; lon <= geographical_area.lr.x; lon += .05) {
		for (double fudge = 0; fudge <= 5; fudge += 0.5) {
			projection->forward(lon, geographical_area.ul.y-fudge, 
					    &(projected.x), &(projected.y));
			updateMinbox(projected.x, projected.y, &projarea);
		}

	}

	// Check bottom of map
	for (double lon = geographical_area.ul.x; lon <= geographical_area.lr.x; lon += .05) {
		projection->forward(lon, geographical_area.lr.y, &(projected.x), &(projected.y));
		updateMinbox(projected.x, projected.y, &projarea);
	}

	// Check left of map
	for (double lat = geographical_area.ul.y; lat >= geographical_area.lr.y; lat -= .05) {
		projection->forward(geographical_area.ul.x, lat, &(projected.x), &(projected.y));
		updateMinbox(projected.x, projected.y, &projarea);
	}

	// Check right of map
	for (double lat = geographical_area.ul.y; lat >= geographical_area.lr.y; lat -= .05) {
		projection->forward(geographical_area.lr.x, lat, &(projected.x), &(projected.y));
		updateMinbox(projected.x, projected.y, &projarea);
	}


	// Check a bunch of points
	for (double lon = geographical_area.ul.x; lon < geographical_area.lr.x; lon += .5) {
		for (double lat = geographical_area.lr.y; lat < geographical_area.ul.y; lat += .5) {
			projection->forward(lon, lat, &(projected.x), &(projected.y));
			updateMinbox(projected.x, projected.y, &projarea);
		}
	}

	return projarea;
}

Area FindRasterExtent(shared_ptr<ProjectedRaster> raster,
		      Area geographical_area)
{
	Area rasterArea;
	double pixel_size = raster->pixel_size;
	shared_ptr<Projection> projection(raster->getProjection());
	Area projArea = FindProjectedExtent(projection, geographical_area, raster->pixel_size);
	int rows = (projArea.lr.x - projArea.ul.x) / pixel_size;
	int cols = (projArea.ul.y - projArea.lr.y) / pixel_size;
	int first_row = (raster->ul_x - projArea.ul.x) / pixel_size;
	int first_col = 0;

	rasterArea.ul.y = first_row;
	rasterArea.ul.x = first_col;
	rasterArea.lr.y = first_row + rows;
	rasterArea.lr.x =  cols;

	return rasterArea;

}
