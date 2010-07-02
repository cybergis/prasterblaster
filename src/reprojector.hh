/* 
   Programmer: David Mattli
*/

#ifndef REPROJECTOR_HH
#define REPROJECTOR_HH

#include <vector>


#include "pRPL/prProcess.h"
#include "pRPL/neighborhood.h"
#include "pRPL/cellSpace.h"
#include "pRPL/layer.h"

//#include "gctp_cpp/projection.h"

#include "projectedraster.hh"
#include "resampler.hh"

using namespace pRPL;
using resampler::resampler_func;

class Reprojector
{
public:
	Reprojector(PRProcess prc, 
		    ProjectedRaster *_input, 
		    ProjectedRaster *_output);
	~Reprojector();
	void reproject();
	void parallelReproject();


private:

	vector<vector<int> > getIndices(int size, int count);
	void reprojectChunk(int firstRow, int numRows);
	double maxx, minx, maxy, miny;
	ProjectedRaster *input;
	ProjectedRaster *output;
	PRProcess prc;
	Layer<unsigned char> input_layer;
	Layer<unsigned char> output_layer;


	resampler_func resampler;
};
void FindMinBox2(double in_ul_x, double in_ul_y,
		 double in_pix_size,
		 int rows, int cols, 
		 Projection *inproj,
		 Projection *outproj,
		 double out_pixsize,
		 double &_ul_x, double &_ul_y, double &_lr_x, double &_lr_y);

void FindMinBox(ProjectedRaster *input, Projection *outproj,
		double out_pixsize,
		double &_ul_x, double &_ul_y, double &_lr_x, double &_lr_y);

ProjectedRaster* GetOutputRaster(ProjectedRaster* input,
				 Projection *out_proj,
				 double out_pixel_size);

#endif // REPROJECTOR_HH
