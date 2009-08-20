

#ifndef REPROJECTOR_HH
#define REPROJECTOR_HH

#include <vector>

#include "gctp_cpp/projection.h"

#include "projectedraster.hh"

class Reprojector
{
public:
  Reprojector(ProjectedRaster *_input, ProjectedRaster *_output);
  ~Reprojector();
  void reproject();
  void parallelReproject(int rank, int numProcs);

private:
  //  void findMinBox();
  double maxx, minx, maxy, miny;
  ProjectedRaster *input;
  ProjectedRaster *output;
  void (*resampler)(void *inraster, double in_x,
		      double in_y, unsigned long in_cols,
		      void *outraster,
		      unsigned long out_x, unsigned long out_y,
		      unsigned long out_cols);
};

void FindMinBox(ProjectedRaster *input, Projection *outproj,
		double out_pixsize,
		double &_ul_x, double &_ul_y, double &_lr_x, double &_lr_y);



#endif // REPROJECTOR_HH
