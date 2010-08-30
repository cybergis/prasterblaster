/* 
   Programmer: David Mattli
*/

#ifndef REPROJECTOR_HH
#define REPROJECTOR_HH

#include <vector>

//#include "gctp_cpp/projection.h"

#include "projectedraster.hh"
#include "resampler.hh"

using resampler::resampler_func;

class Reprojector
{
public:
	Reprojector(ProjectedRaster *_input, 
		    ProjectedRaster *_output);
	~Reprojector();
	void reproject();
	void parallelReproject();


private:
	void reprojectChunk(int firstRow, int numRows);
	int numprocs, rank;
	double maxx, minx, maxy, miny;
	ProjectedRaster *input;
	ProjectedRaster *output;

	resampler_func resampler;
};
Area FindMinBox2(double in_ul_x, double in_ul_y,
		 double in_pix_size,
		 int rows, int cols, 
		 Projection *inproj,
		 Projection *outproj,
		 double out_pixsize);


void FindMinBox(ProjectedRaster *input, Projection *outproj,
		double out_pixsize,
		double &_ul_x, double &_ul_y, double &_lr_x, double &_lr_y);

ProjectedRaster* GetOutputRaster(ProjectedRaster* input,
				 Projection *out_proj,
				 double out_pixel_size);

#endif // REPROJECTOR_HH
