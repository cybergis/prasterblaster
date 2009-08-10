
#include <cstdio>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/coordinate.h"

#include "reprojector.hh"
#include "resampler.hh"

Reprojector::Reprojector(ProjectedRaster *_input, ProjectedRaster *_output) :
  input(_input), output(_output)
{
  maxx = maxy = 0;
  minx = miny = 1e+37;

  if (input->getBitCount() == 8) {
    resampler = &nearest_neighbor<unsigned char>;
  }

  return;
};

Reprojector::~Reprojector()
{
  
  
  return;
}
void parallelReproject(int rank, int numProcs)
{
  
  return;
}
void Reprojector::reproject()
{
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Assume output memory is allocated!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Transformer t;
  Coordinate temp, temp2;
  double in_ulx, in_uly, out_ulx, out_uly;
  double in_pixsize, out_pixsize;
  long in_rows, in_cols, out_rows, out_cols;


  t.setInput(*output->getProjection());
  t.setOutput(*input->getProjection());
  
  in_rows = input->getRowCount();
  in_cols = input->getColCount();
  out_rows = output->getRowCount();
  out_cols = output->getColCount();
  in_pixsize = input->getPixelSize();
  out_pixsize = output->getPixelSize();
  in_ulx = input->ul_x;
  in_uly = input->ul_y;
  out_ulx = output->ul_x;
  out_uly = output->ul_y;
  

  temp.units = METER;
  for (int y = 0; y < out_rows; ++y) {
    if (y % 100 == 0)
      printf("On col %d of %ld\n", y, out_rows);
    for (int x = 0; x < out_cols; ++x) {
      temp.x = ((double)x * out_pixsize) + out_ulx;
      temp.y = ((double)y * out_pixsize) - out_uly;
      t.transform(&temp);
      output->getProjection()->inverse(temp.x, temp.y);
	output->getProjection()->forward(output->getProjection()->lon(),
					 output->getProjection()->lat(), &(temp2.x), &(temp2.y));
	temp.x = ((double)x * out_pixsize) + out_ulx;
	temp.y = ((double)y * out_pixsize) - out_uly;
	
	t.transform(&temp);
	
	// temp now contains coords to input projection
	temp.x -= in_ulx;
	temp.y += in_uly;
	temp.x /= in_pixsize;
	temp.y /= out_pixsize;
	//      printf("temp.x: %lld, temp.y: %lld, x: %d, y: %d\n", 
	//	     (long long)temp.x, (long long)temp.y, x, y);
	// temp is now scaled to input raster coords, now resample!
	resampler(input->data, temp.x, temp.y, in_cols, 
		  output->data, x, y, out_cols);

    }
  }

  return;
}


// Minbox finds number of rows and columns in output and upper-left
// corner of output
void FindMinBox(ProjectedRaster *input, Projection *outproj, double out_pixsize,
		double &_ul_x, double &_ul_y, double &_lr_x, double &_lr_y)
{
  double ul_x, ul_y;
  double pixsize; // in METER!
  Coordinate temp;
  Transformer t;
  int cols, rows;
  double maxx, minx, maxy, miny;
  double ul_lon, ul_lat, lr_lon, lr_lat;

  maxx = maxy = 0;
  minx = miny = 1e+37;

  ul_x = input->ul_x;
  ul_y = input->ul_y;
  cols = input->cols;
  rows = input->rows;
  pixsize = input->getPixelSize();
  
  t.setInput(*input->getProjection());
  t.setOutput(*outproj);

  // Find geographic corners of input
  input->getProjection()->inverse(ul_x, ul_y, &ul_lon, &ul_lat);
  input->getProjection()->inverse(ul_x+(cols*pixsize), ul_y-(rows*pixsize), &lr_lon,  &lr_lat);
  double delta_east = (lr_lon-ul_lon)/cols, delta_north = (ul_lat-lr_lat)/rows;

  // Calculate minbox
  temp.x = ul_x;
  temp.y = ul_y;
  temp.units = METER;

  // Check top of map
  for (int x = 0; x < cols; ++x) {
    temp.x = (double)x*pixsize + ul_x;
    temp.y = ul_y;
    t.transform(&temp);
    if (temp.x < minx) minx = temp.x;
    if (temp.x > maxx) maxx = temp.x;
    if (temp.y < miny) miny = temp.y;
    if (temp.y > maxy) maxy = temp.y;
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
    temp.x = (double)x*pixsize + ul_x;
    temp.y = (double)-rows*pixsize + ul_y;
    t.transform(&temp);
    if (temp.x < minx) minx = temp.x;
    if (temp.x > maxx) maxx = temp.x;
    if (temp.y < miny) miny = temp.y;
    if (temp.y > maxy) maxy = temp.y;
    temp.x = ul_lon + (x * delta_east);
    temp.y = ul_lat + (rows*delta_north);
    outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
    if (temp.x < minx) minx = temp.x;
    if (temp.x > maxx) maxx = temp.x;
    if (temp.y < miny) miny = temp.y;
    if (temp.y > maxy) maxy = temp.y;

  }
 
  // Check Left side
  for (int y = 0; y < rows; ++y) {
    temp.x = ul_x;
    temp.y = (double)-y*(pixsize+1) + ul_y;
    t.transform(&temp);
    if (temp.x < minx) minx = temp.x;
    if (temp.x > maxx) maxx = temp.x;
    if (temp.y < miny) miny = temp.y;
    if (temp.y > maxy) maxy = temp.y;
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
    temp.x = (double)cols*(1+pixsize) + ul_x;
    temp.y = (double)-y*pixsize + ul_y;
    t.transform(&temp);
    if (temp.x < minx) minx = temp.x;
    if (temp.x > maxx) maxx = temp.x;
    if (temp.y < miny) miny = temp.y;
    if (temp.y > maxy) maxy = temp.y;
    temp.x = ul_lon + (cols * delta_east);
    temp.y = ul_lat + (y * delta_north);
    outproj->forward(temp.x, temp.y, &(temp.x), &(temp.y));
    if (temp.x < minx) minx = temp.x;
    if (temp.x > maxx) maxx = temp.x;
    if (temp.y < miny) miny = temp.y;
    if (temp.y > maxy) maxy = temp.y;
  }


  printf("Min x: %e, Max x: %e, Max y: %e, Min y: %e\n",
	 minx, maxx, maxy, miny);
  printf("Output Raster Dim: %lld cols, %lld rows\n", 
	 (long long)((maxx-minx)/out_pixsize),
	 (long long)((maxy-miny)/out_pixsize));

  // Set outputs
  _ul_x = minx;
  _ul_y = maxy;
  _lr_x = maxx;
  _lr_y = miny;
  

  return;
}

