/*
  Programmer: David Mattli
*/

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cctype>

#include <mpi.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mollweide.h"
#include "projectedraster.hh"
#include "reprojector.hh"

#include <gdal_priv.h>

#include "driver.hh"

const char *usage = "usage: prasterblaster <input raster path> <output raster path> " 
  "<output raster projection> \n"; 


int main(int argc, char *argv[]) 
{
 	long int rows, cols;

	ProjectedRaster *out;
	int rank;

        rows = cols = 0;

	MPI_Init(&argc, &argv);
	if (argc < 4) {
		if (rank == 0) {
			printf("%s\n", usage);
		}
		MPI_Abort(MPI_COMM_WORLD, -1);
		return 0;
        }

	driver(argv[1], argv[2], argv[3]);
	MPI_Finalize();
	
	return 0;
}
