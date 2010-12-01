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

char *usage = "usage: prasterblaster <input raster path> <output raster path> " 
  "<out raster description path(xml)>\n";


int main(int argc, char *argv[]) 
{
 	long int rows, cols;

	ProjectedRaster *out;
	Reprojector *re = 0;
	int rank;

        rows = cols = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank); 

	if (argc < 4) {
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


	if (rank == 0) {
		out = new ProjectedRaster(argv[2],
					  argv[3]);

		out = new ProjectedRaster(argv[2], 
					  &in,
					  
		delete out;
		out = 0;
		MPI_Barrier(MPI_COMM_WORLD);
		printf("Output raster created and nodes synced!\n");
	} else {
		MPI_Barrier(MPI_COMM_WORLD);
	}


	// Now we re-open the output raster on each node.
	out = new ProjectedRaster(argv[2]);
	if (out == 0) {
		fprintf(stderr, "Output allocation failed, something is very wrong!\n");
		MPI_Finalize();
		return 1;

	}
	
	if (!in.isReady()) {
		fprintf(stderr, "Error opening input raster, not ready!\n");
		MPI_Finalize();
		return 1;

	}

	if (!out->isReady()) {
		fprintf(stderr, "Error opening output raster, not ready!\n");
		MPI_Finalize();
		return 1;

	}

	re = new Reprojector(&in, out);
	re->parallelReproject();

	// Cleanup
	delete re;
	delete out;

	MPI_Finalize();
	
	return 0;
}
