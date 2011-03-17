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
		return value;
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
Reprojector::Reprojector(shared_ptr<ProjectedRaster> _input, shared_ptr<ProjectedRaster> _output,
			 int processCount, int rank) 
        :
        input(_input), output(_output), numprocs(processCount), rank(rank)
{
	numchunks = numprocs * 2;
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
	int chunk_count = numprocs * 2;
	vector<ChunkExtent> chunks = output->getChunks(chunk_count);
	vector<int> assignments = getAssignments(chunk_count);

	for (unsigned int i = 0; i < assignments.size(); ++i) {
		reprojectChunk(chunks[assignments[i]]);
	}


  /*
	int chunk_count = numprocs * 2;
	int process_count = numprocs;
	std::vector<ChunkExtent> ce = getChunkExtents(output->getRowCount(), 
						      chunk_count, process_count);

	for (int i=0; i < ce.size(); ++i) {
		if (ce[i].processAssignment() == rank) {
			reprojectChunk(ce[i].firstIndex(), (ce[i].lastIndex())-(ce[i].firstIndex()));
		}
	}
	

        return;
  */
}

vector<int> Reprojector::getAssignments(int chunkCount)
{
	vector<int> assign;
	int normal_count = chunkCount / numprocs;
	int last_count = normal_count + (chunkCount % numprocs);
	int count = normal_count;

	// TODO: Handle case where numprocs < chunkCount better
	
	if (rank == (numprocs -1)) {
		count = last_count;
	}
	
	for (int i = normal_count * rank; i < count; ++i) {
		assign.push_back(i);
	}

	return assign;
}

void Reprojector::reprojectChunk(ChunkExtent chunk)
{
        shared_ptr<Projection> outproj, inproj;
        Coordinate temp1, temp2;
        Coordinate in_ul, in_lr, out_ul;
        double in_pixsize, out_pixsize; 
	long in_rows, in_cols, out_rows, out_cols;
	std::vector<char> inraster, outraster;
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
	
	input_proj_area = FindProjectedExtent(inproj,
					      chunk.getGeographicalMinbox(),
					      in_pixsize);
					      
	
        in_ul.x = input->ul_x;
	in_ul.y = input_proj_area.ul.y;
	in_lr.x = input_proj_area.lr.x;
	in_lr.y = input_proj_area.lr.y;
	
        if (in_ul.y > input->ul_y)
                in_ul.y = input->ul_y;
        in_lr.x = input->ul_x + (in_pixsize * in_cols);
	
        long in_first_row = (long)((input->ul_y - in_ul.y) / in_pixsize);
        in_rows = (long)ceil((in_ul.y - in_lr.y) / in_pixsize);


        // Setup raster vectors
        size_t s = out_rows;
        s *= out_cols;
        s *= output->bitsPerPixel()/8;
        outraster.resize(s);
        s = in_rows;
        s *= in_cols;
        s *= input->bitsPerPixel()/8;
        inraster.resize(s);
	
	printf("Reading input...\n");
        // Read input file
        if (input->readRaster(in_first_row, in_rows, &(inraster[0]))) {
                //              printf("Read %d rows\n", numRows);
        } else {
		fprintf(stderr, "Error Reading input!\n");
		fflush(stderr);
		return;
        }
	
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
			pixelArea = rt.Transform(temp1);   /// !!! NO No, this is bad. Makes no sense.
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


			for (int i=0, a=0; i < pixel_height; ++i) {
				for (int j=0; j < pixel_width; ++j) {
					index_array[++a] = (j+(long)temp1.x) 
						+ ((i+(long)temp2.y) *in_cols);
				}
			}
			temp1.x = floor(temp1.x);
			temp1.y = floor(temp1.y);
			unsigned int index = 0;
			try {
				index = (unsigned int)temp1.x;
				index += ((unsigned int)temp1.y) * in_cols;

				outraster.at(chunk_x+(chunk_y*out_cols)) = inraster.at(index);
				++total;
			} catch (out_of_range &e) {
				printf("------\n");
				printf("Attemped X: %f Y: %f, Raster dimensions X: %ld, Y: %ld\n",
				       temp1.x, temp1.y,
				       in_cols, in_rows);
				outraster.at(chunk_x+(chunk_y*out_cols)) = 255;
				++count;
				printf("INDEX: %d size: %u\n", index, inraster.size());

			}
		
		}
	}
	printf("%d total, %d pixels off!\n", total, count);
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

Area FindGeographicalExtent(shared_ptr<Projection> projection, 
			    Coordinate ul_point,
			    int row_count,
			    int col_count,
			    double pixel_size)
{
	Area geoarea;
	double lon, lat;
	double x, y;
	
	geoarea.ul.x = DBL_MAX;
	geoarea.ul.y = -DBL_MAX;
	geoarea.lr.x = -DBL_MAX;
	geoarea.lr.y = DBL_MAX;
	lon = lat = x = y = 0.0;

	for (int row = 0; row < row_count; ++row) {
		for (int col = 0; col < col_count; ++col) {
			x = ul_point.x + (pixel_size * col);
			y = ul_point.y - (pixel_size * row);
			projection->inverse(x, y, &lon, &lat);
			
			updateMinbox(lon, lat, &geoarea);
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

	// Check top of map
	for (double lon = geographical_area.ul.x; lon <= geographical_area.lr.x; lon+=.05) {
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
	for (double lon = -180.0; lon < 180.0; lon += .5) {
		for (double lat = -90.0; lat < 90.0; lat += .5) {
			projection->forward(lon, lat, &(projected.x), &(projected.y));
			updateMinbox(projected.x, projected.y, &projarea);
		}
	}
	

	return projarea;
}
      
			    

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

