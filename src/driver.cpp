/*
  Programmer: David Mattli
*/

#include <cstdio>
#include <cstring>


#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mollweide.h"
#include "projectedraster.hh"
#include "reprojector.hh"
#include "rasterreader.hh"

#include <gdal_priv.h>


double params[15] =  { 6370997.000000, 
		       0, 0, 0, 0, 
		       0, 0, 0, 0, 0, 
		       0, 0, 0, 0, 0};



int main(int argc, char *argv[]) 
{
	double ul_x, ul_y, lr_x, lr_y;
 	long int rows, cols;

	vector<unsigned char> *dat = 0;
	PRProcess prc(MPI_COMM_WORLD);
	ProjectedRaster *out;
	Reprojector *re = 0;

        rows = cols = 0;

	prc.init(argc, argv);


	if (argc < 3) {
                if (prc.isMaster())
                        printf("usage: prasterblaster <input raster path> <output raster path>\n");
                prc.abort();
                return(1);
        }
	

        ProjectedRaster in(argv[1]);

        if (in.isReady() == true) {
                if(prc.isMaster())
                        printf("Input raster opened.\n");
        } else {
                
                prc.abort();
                return 1;
        }

	Projection *outproj;
	outproj = new Mollweide(params, METER, in.getDatum());

	if (prc.isMaster()) {
		out = new ProjectedRaster(argv[2],
					  &in,
					  outproj,
					  in.getPixelType(),
					  in.getPixelSize());

		delete out;
		out = 0;
		prc.sync();
		printf("synced!\n");
	} else {
		prc.sync();
	}

	out = new ProjectedRaster(argv[2]);

	if (out != 0 && !(in.isReady() && out->isReady())) {
		printf("Error in opening rasters\n");
		prc.finalize();
		return 1;
	}
	

	re = new Reprojector(prc, &in, out);
	re->parallelReproject();

	// Cleanup
	delete re;
	delete out;
	delete outproj;	

	prc.finalize();
	
	return 0;
}
