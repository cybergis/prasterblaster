/*
  Programmer: David Mattli
*/

#include <boost/mpi.hpp>
#include <mpi.h>


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

double params[15] =  { 6370997.000000, 
		       0, 0, 0, 0, 
		       0, 0, 0, 0, 0, 
		       0, 0, 0, 0, 0};



int main(int argc, char *argv[]) 
{
	double ul_x, ul_y, lr_x, lr_y;
 	long int rows, cols;
	mpi::environment env(argc, argv);

	// mpi::communicator world;

	/*
	rows = cols = 0;
	ProjectedRaster in("/home/dmattli/Desktop/mmr/veg_geographic_1deg.img");
	if (in.isReady() == true) {
		printf("Image read!\n");
	}

	Projection *outproj;
	outproj = new Hammer(params, METER, in.getDatum());

	FindMinBox(&in, outproj, in.bitsPerPixel(), ul_x, ul_y, lr_x, lr_y);
	FindMinBox(&in, outproj, in.bitsPerPixel(), ul_x, ul_y, lr_x, lr_y);
	rows = (ul_y-lr_y) / in.getPixelSize();
	cols = (lr_x-ul_x) / in.getPixelSize();

	ProjectedRaster out("/home/dmattli/Desktop/output.tif", 
			    in.getRowCount(), in.getColCount(), 
			    in.getPixelType(), in.getPixelSize(), 
			    in.bandCount(), outproj, ul_x, ul_y);

	//x	Reprojector rp(&in, &out); 
	//	rp.reproject();

	if (world.rank() == 0) {

	}


	*/
	return 0;
}
