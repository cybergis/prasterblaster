/*
  Programmer: David Mattli
*/

#include <cstdio>
#include <cstring>

#include <mpi.h>

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

char *usage = "usage: prasterblaster <input raster path> <output raster path>\n";



int main(int argc, char *argv[]) 
{
	double ul_x, ul_y, lr_x, lr_y;
 	long int rows, cols;

	vector<unsigned char> *dat = 0;
	ProjectedRaster *out;
	Reprojector *re = 0;
	int rank;

        rows = cols = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank); 

	if (argc < 3) {
		if (rank == 0) {
			printf(usage);
		}
		MPI_Abort(MPI_COMM_WORLD, -1);
		return 0;
        }
	
	
        ProjectedRaster in(argv[1]);

        if (in.isReady() == true) {
                if(rank == 0)
                        printf("Input raster opened.\n");
        } else {
                
		MPI_Abort(MPI_COMM_WORLD, -1);
                return 1;
        }

	Projection *outproj;
	outproj = new Mollweide(params, METER, in.getDatum());

	if (rank == 0) {
		out = new ProjectedRaster(argv[2],
					  &in,
					  outproj,
					  in.getPixelType(),
					  in.getPixelSize());

		delete out;
		out = 0;
		MPI_Barrier(MPI_COMM_WORLD);
		printf("synced!\n");
	} else {
		MPI_Barrier(MPI_COMM_WORLD);
	}

	out = new ProjectedRaster(argv[2]);
	if (out == 0) {
		fprintf(stderr, "Error reopening output raster\n");
		MPI_Finalize();
		return 1;

	}
	
	if (!in.isReady()) {
		fprintf(stderr, "Error opening input raster!\n");
		MPI_Finalize();
		return 1;

	}

	if (!out->isReady()) {
		fprintf(stderr, "Error opening output raster!\n");
		MPI_Finalize();
		return 1;

	}

	re = new Reprojector(&in, out);
	re->parallelReproject();

	// Cleanup
	delete re;
	delete out;
	delete outproj;	

	MPI_Finalize();
	
	return 0;
}
