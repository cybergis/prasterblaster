
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>


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


int main(int argc, char *argv[])
{

  double ul_x, ul_y, lr_x, lr_y;
  long int rows, cols;
  rows = cols = 0;

    ProjectedRaster *in = RasterReader::readImgRaster("/home/dmattli/Desktop/holdnorm_geographic_30min");
  //  ProjectedRaster *in = RasterReader::readRaster("/home/dmattli/Desktop/input.tif");
    if (in == 0) {
      printf("Oh no!\n");
    exit(-1);
    } else {
      printf("Image read!\n");
    }

  Projection *outproj;
  outproj = new Hammer(params, METER, (ProjDatum)19);
  
  FindMinBox(in, outproj, in->getPixelSize(), ul_x, ul_y, lr_x, lr_y);
  cols = (int)((lr_x-ul_x)/in->getPixelSize());
  rows = (int)((ul_y-lr_y)/in->getPixelSize());
  
  ProjectedRaster out(rows, cols, in->getBitCount());
  out.setProjection(HAMMER);
  out.setDatum((ProjDatum)19);
  out.setGctpParams(params);
  out.setUL(ul_x, ul_y);
  out.setPixelSize(in->getPixelSize());
  out.setUnit(METER);

  Reprojector rp(in, &out);
  rp.reproject();

  //  RasterReader::readRaster("/home/dmattli/Desktop/input.tif");
  RasterReader::writeRaster(string("/home/dmattli/Desktop/puppies.tif"),
			    &out);
  delete outproj;
  delete in;
  return 0;
}
