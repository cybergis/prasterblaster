
#include <boost/mpi.hpp>


#include <QImage>
#include <QString>

#include <cstdio>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/hammer.h"
#include "projectedraster.hh"
#include "reprojector.hh"
#include "rasterreader.hh"

#include <gdal_priv.h>

namespace mpi = boost::mpi;

double params[15] =  { 6370997.000000, 0, 00000000.000000, 00000000.000000, 0.000000, 
		       0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
		       0.000000, 0.000000, 0.000000, 0.000000, 0.000000};



int main(int argc, char *argv[]){


	double ul_x, ul_y, lr_x, lr_y;
	long int rows, cols;
	mpi::environment env(argc, argv);
	mpi::communicator world;

	rows = cols = 0;
	ProjectedRaster in("/home/dmattli/Desktop/holdnorm_geographic_30min");
	if (in.isReady() == true) {
		printf("Image read!\n");
	}
	Projection *outproj;
	outproj = new Hammer(params, METER, (ProjDatum)19);
  
	FindMinBox(&in, outproj, in.getPixelSize(), ul_x, ul_y, lr_x, lr_y);
	cols = (int)((lr_x-ul_x)/in.getPixelSize());
	rows = (int)((ul_y-lr_y)/in.getPixelSize());

	ProjectedRaster out(rows, cols, in.getBitCount(), outproj, ul_x, ul_y);
	out.setUnit(METER);

	Reprojector rp(&in, &out);
	//  rp.reproject();
	rp.parallelReproject(world.rank(), world.size());
	if (world.rank() == 0)
		out.write("/home/dmattli/Desktop/output.tif");

	return 0;
}
