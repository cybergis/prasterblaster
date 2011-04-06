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
#include <vector>

#include <boost/shared_ptr.hpp>

#include <mpi.h>
#include <gdal.h>
#include <gdal_priv.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/coordinate.h"

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

	if (fabs(temp1.x - temp2.x) > 0.0001) {
		// Overlap detected!
		throw std::string("Point Overlaps!");
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
	Coordinate chunk_ul;
	Area geominbox, projminbox;

	if (dest_raster->getRowCount() <= 0 || process_count <= 0) {
		return chunks;
	}

	for (int i=0; i < process_count-1; ++i) {
		chunk_ul.x = source_raster->ul_x;
		chunk_ul.y = source_raster->ul_y - i * chunk_size * dest_raster->pixel_size;
		geominbox = FindGeographicalExtent(dest_raster->getProjection(),
						   chunk_ul,
						   chunk_size,
						   dest_raster->getColCount(),
						   dest_raster->pixel_size);
		clampGeoCoordinate(&geominbox.ul);
		clampGeoCoordinate(&geominbox.lr);
		projminbox = FindProjectedExtent(source_raster->getProjection(),
						 geominbox,
						 source_raster->pixel_size);

		chunks.push_back(ChunkExtent(i*chunk_size,
					     ((i+1)*chunk_size-1),
					     geominbox,
					     projminbox));

		
	}

	// Last chunk includes any leftover
	long last_row = source_raster->getRowCount() - 1;
	chunk_ul.x = source_raster->ul_x;
	chunk_ul.y = source_raster->ul_y - ((process_count - 1) 
					    * chunk_size * source_raster->pixel_size);
	long last_chunk_size = last_row - (process_count - 1) + 1;
	geominbox = FindGeographicalExtent(source_raster->getProjection(),
					   chunk_ul,
					   last_chunk_size,
					   source_raster->getColCount(),
					   source_raster->pixel_size);
	clampGeoCoordinate(&geominbox.ul);
	clampGeoCoordinate(&geominbox.lr);
	projminbox = FindProjectedExtent(source_raster->getProjection(),
					 geominbox,
					 dest_raster->pixel_size);
	

	chunks.push_back(ChunkExtent((process_count-1)*chunk_size, last_row, geominbox,
				 projminbox));

	return chunks;

	
}
	

Reprojector::Reprojector(shared_ptr<ProjectedRaster> _input, shared_ptr<ProjectedRaster> _output,
			 int processCount, int rank) 
        :
        input(_input), output(_output), numprocs(processCount), rank(rank)
{
	numchunks = numprocs * 100;
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
	printf("Finding chunks\n");

	printf("Getting assignments\n");


}

void Reprojector::reprojectChunk(ChunkExtent chunk)
{

	printf("Reprojecting chunk %d to %d...", chunk.firstIndex(), chunk.lastIndex());
	fflush(stdout);
        shared_ptr<Projection> outproj, inproj;
        Coordinate temp1, temp2;
        Coordinate in_ul, in_lr, out_ul;
        double in_pixsize, out_pixsize; 
	long in_rows, in_cols, out_rows, out_cols;
	std::vector<char> inraster, outraster;
	set<long> overflow_rows;
	Area input_proj_area;
	int firstRow = chunk.firstIndex();
	int numRows = chunk.lastIndex() - chunk.firstIndex() + 1;

        if (firstRow + numRows > output->getRowCount()) {
                fprintf(stderr, "Invalid chunk range... %d is > %d\n",
                        firstRow + numRows,
                        output->getRowCount());
                fflush(stderr);
                return;
        }

        out_pixsize = output->getPixelSize();
        in_pixsize = input->getPixelSize();
        outproj = shared_ptr<Projection>(output->getProjection());
        inproj = shared_ptr<Projection>(input->getProjection());

        out_rows = numRows;
        out_cols = output->getColCount();
        out_ul.x = output->ul_x;
        out_ul.y = output->ul_y - firstRow * out_pixsize;
        in_cols = input->getColCount();
	
	Area in_raster_area = FindRasterArea(input, output, chunk.firstIndex(),
					     chunk.lastIndex());
	in_ul = input->getProjectedCoordinate(0, in_raster_area.ul.y);
	in_lr = input->getProjectedCoordinate(in_raster_area.lr.x, in_raster_area.lr.y);
	long in_first_row = in_raster_area.ul.y;
	in_rows = in_raster_area.ul.y - in_lr.y + 1;
	

	printf("First row: %d, count %d\n", in_first_row, in_rows);

        // Setup raster vectors
        size_t s = out_rows;
        s *= out_cols;
        s *= output->bitsPerPixel()/8;
	printf("Allocating %u MB of output memory, %d rows\n", s/1024/1024, out_rows);
        outraster.resize(s);
        s = in_rows + 1;
        s *= in_cols;
        s *= input->bitsPerPixel()/8;
	printf("Allocating %u MB of input memory, %d rows\n", s/1024/1024, in_rows);
        inraster.resize(s);

	printf("Reading input...");
	fflush(stdout);
        // Read input file
	try {
		input->readRaster(in_first_row, in_rows, &(inraster[0]));
	} catch (std::exception) {
		fprintf(stderr, "Error Reading input!\n");
		fflush(stderr);
		return;
        }
	printf("done\n");
	Area pixelArea;
	int count = 0;
	int total = 0;
	RasterCoordTransformer rt(output, input);        

	for (int chunk_y = 0; chunk_y < out_rows; ++chunk_y)  {
		for (int chunk_x = 0; chunk_x < out_cols; ++chunk_x) {

			// Determine location of equivalent input pixel
			outraster.at(chunk_y*out_cols) = 127; // REMOVE THIS

			temp1.x = chunk_x; 
			temp1.y = chunk_y + firstRow;
			pixelArea = rt.Transform(temp1);   /// Now makes sense.
			temp1 = pixelArea.ul;
			temp2 = pixelArea.lr;

			//! TODO: Check for overlap!
			// temp1&2 are now scaled to input raster coords, now resample!
			// But does the rectangle defined by temp1 and temp2 actually
			// contain any points? If not use nearest-neighbor...
			long center_index = (long)((temp1.x + temp2.x)/2);
			center_index += (long)((temp1.y + temp2.y)/2) * in_cols;
			long pixel_width = (long)fabs((temp2.x - temp1.x)/in_pixsize);
			long pixel_height = (long)(fabs(temp1.y - temp2.y)/in_pixsize);

			if (pixel_width == 0) {
				pixel_width = 1;
			}
			if (pixel_height == 0) {
				pixel_height = 1;
			}
			std::vector<long> index_array(pixel_width*pixel_height, 0L);

/*
			for (int i=0, a=0; i < pixel_height; ++i) {
				for (int j=0; j < pixel_width; ++j, ++a) {
				  index_array.at(a) = (j+(long)temp1.x) 
						+ ((i+(long)temp2.y) *in_cols);
				}
			}
*/
			temp1.x = floor(temp1.x);
			temp1.y = floor(temp1.y);
			unsigned int index = 0;
			try {
				index = (unsigned int)temp1.x;
				index += ((unsigned int)temp1.y) * in_cols;

				outraster.at(chunk_x+(chunk_y*out_cols)) = inraster.at(index);
				++total;
			} catch (out_of_range &e) {
				overflow_rows.insert(temp1.y);
				outraster.at(chunk_x+(chunk_y*out_cols)) = 255;
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
	output->writeRaster(firstRow, numRows, &(outraster[0]));
	
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
			if (fabs(temp.x - lon) > 0.0001) {
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
			if (fabs(temp.x - lon) > 0.0001) {
				continue;
			} else {
				updateMinbox(lon, lat, &geoarea);
			}
		}
	}

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
		projection->forward(lon, geographical_area.ul.y, &(projected.x), &(projected.y));
		updateMinbox(projected.x, projected.y, &projarea);
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

	return rasterArea;}

      
			    

Area FindMinBox(double in_ul_x, double in_ul_y,
		 double in_pix_size,
		 int rows, int cols, 
		 Projection *inproj,
		 Projection *outproj,
		 double out_pixsize)
{
	double ul_x, ul_y;
	double pixsize; // in METER!
	Coordinate temp;
	Transformer t;
	double maxx, minx, maxy, miny;
	double ul_lon, ul_lat, lr_lon, lr_lat;
	Area area;

	maxx = maxy = 0;
	minx = miny = 1e+37;

	ul_x = in_ul_x;
	ul_y = in_ul_y;
	pixsize = in_pix_size;
	pixsize = out_pixsize;
	
  
	t.setInput(*(inproj->copy()));
	t.setOutput(*(outproj->copy()));

	// Find geographic corners of input
	inproj->inverse(ul_x, ul_y, &ul_lon, &ul_lat);
	inproj->inverse(ul_x+(cols*pixsize), 
					ul_y-(rows*pixsize), 
					&lr_lon, 
					&lr_lat);
	double delta_east = (lr_lon-ul_lon)/cols, 
		delta_north = (ul_lat-lr_lat)/rows;

	// Calculate minbox
	temp.x = ul_x;
	temp.y = ul_y;
	temp.units = METER;

	// Check top of map
	for (int x = 0; x < cols; ++x) {
	  /*		temp.x = (double)x*pixsize + ul_x;
		temp.y = ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y; */
		temp.x = ul_lon + (x*delta_east);
		temp.y = ul_lat; 
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
	}

	// Check bottom of map
	for (int x = 0; x < cols; ++x) {
	  /*		temp.x = (double)x*pixsize + ul_x;
		temp.y = (double)-rows*pixsize + ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;*/
		temp.x = ul_lon + (x * delta_east);
		temp.y = ul_lat - (rows*delta_north); 
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;

	}
 
	// Check Left side
	for (int y = 0; y < rows; ++y) {
	  /*		temp.x = ul_x;
		temp.y = (double)-y*(pixsize+1) + ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y; */
		temp.x = ul_lon;
		temp.y = ul_lat - (y*delta_north);
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
	}

	// Check right side
	for (int y = 0; y < rows; ++y) {
	  /*		temp.x = (double)cols*(1+pixsize) + ul_x;
		temp.y = (double)-y*pixsize + ul_y;
		t.transform(&temp);
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y; */
		temp.x = ul_lon + (cols * delta_east); 
		temp.y = ul_lat - (y * delta_north);
		outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
		if (temp.x < minx) minx = temp.x;
		if (temp.x > maxx) maxx = temp.x;
		if (temp.y < miny) miny = temp.y;
		if (temp.y > maxy) maxy = temp.y;
	}

	// Set outputs
	area.ul.x = minx;
	area.ul.y = maxx;
	area.lr.x = maxx;
	area.lr.y = miny;
	


 	outproj->inverse(minx, maxy, &temp.x, &temp.y);
	outproj->inverse(maxx, miny);
//	printf("MINBOX ul %f %f lr %f %f or ul(%f %f) lr(%f %f)\n", minx, maxy, maxx, miny,
//	       temp.x, temp.y, outproj->lon(), outproj->lat());
  

	return area;
}

