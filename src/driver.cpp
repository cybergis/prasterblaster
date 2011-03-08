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
#include "rasterreader.hh"

#include <gdal_priv.h>


double params[15] =  { 6370997.000000, 
		       0, 0, 0, 0, 
		       0, 0, 0, 0, 0, 
		       0, 0, 0, 0, 0};

const char *usage = "usage: prasterblaster <input raster path> <output raster path> " 
  "<output raster projection> \n"; 


int driver(string input_raster, string output_filename, string output_projection)
{
	int rank = 0;
	shared_ptr<ProjectedRaster> in, out;
	shared_ptr<Projection> in_proj, out_proj;
	shared_ptr<Reprojector> re;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

	// Open input raster and check for errors
        in = shared_ptr<ProjectedRaster>(new ProjectedRaster(input_raster));

        if (in->isReady() == true) {
                if (rank == 0)
                        printf("Input raster opened.\n");
        } else {
                fprintf(stderr, "Error opening input raster\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
                return 1;
        }

	// If we are rank 0, create the output projection and raster
	if (rank == 0) {

		// Downcase output projection name
		std::transform(output_projection.begin(),
			       output_projection.end(),
			       output_projection.begin(),
			       (int(*) (int))std::tolower);

		in_proj = shared_ptr<Projection>(in->getProjection());

		if (output_projection == "hammer") {
			out_proj = shared_ptr<Projection>(Transformer::convertProjection(HAMMER));
		} else if (output_projection == "mollweide") {
			out_proj = shared_ptr<Projection>(Transformer::convertProjection(MOLL));
		} else if (output_projection == "sinosoidal") {
			return 1;
		} else {
			// fail...
			return 1;
		}
		out_proj->setUnits(in_proj->units());
		out_proj->setDatum(in_proj->datum());
		out_proj->setParams(in_proj->params());


		bool result  = ProjectedRaster::CreateRaster(output_filename,
							     in,
							     shared_ptr<Projection>(out_proj->copy()),
							     in->type ,
							     in->pixel_size);
		
		if (result == false) {
			fprintf(stderr, "Failed to create output raster!\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		printf("Output raster created and nodes synced!\n");
	} else {
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	
	// Now we re-open the output raster on each node.
	out = shared_ptr<ProjectedRaster>(new ProjectedRaster(output_filename));
	if (out == 0) {
		fprintf(stderr, "Output allocation failed, something is very wrong!\n");
		MPI_Finalize();
		return 1;

	}
	
	if (!in->isReady()) {
		fprintf(stderr, "Error opening input raster, not ready!\n");
		MPI_Finalize();
		return 1;

	}

	if (!out->isReady()) {
		fprintf(stderr, "Error opening output raster, not ready!\n");
		MPI_Finalize();
		return 1;

	}

	re = shared_ptr<Reprojector>(new Reprojector(in, out, 1, 0));
	re->parallelReproject();

	// Cleanup
	

	return 0;
}

int main(int argc, char *argv[]) 
{
 	long int rows, cols;

	ProjectedRaster *out;
	Reprojector *re = 0;
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
