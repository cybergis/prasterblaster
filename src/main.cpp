/*!
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE 
 *
 * This software is in the public domain, furnished "as is", without
 * technical support, and with no warranty, express or implied, as to
 * its usefulness for any purpose.
 *
 * @section DESCRIPTION
 *
 * This file implements a main function of the program, calling
 * MPI_Init, MPI_Finalize and running the driver function.
 *
 */

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include <mpi.h>
#include <getopt.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mollweide.h"
#include "projectedraster.hh"
#include "reprojector.hh"

#include <gdal_priv.h>

#include "driver.hh"

using std::string;

// getopt variables
extern char *optarg;
extern int optind, opterr, optopt;
int analyze_partitions = 0;

const char *usage = "usage: prasterblaster [-n <partition count>] "
	"-p <output projection srs> [-f <fill value>] "
	"<input raster path> <output raster path> \n";

struct option longopts[] = {
	{"analyze-partitions", no_argument, &analyze_partitions, 1},
	{"output-projection", required_argument, NULL, 'p'},
	{"partition-count", required_argument, NULL, 'n'},
	{"resampler", required_argument, NULL, 'r'},
	{"fill-value", required_argument, NULL, 'f'},
	{0, 0, 0, 0}

};


int main(int argc, char *argv[]) 
{
 	long int rows, cols;

	ProjectedRaster *out;
	int rank, c;
	string output_srs, resampler, fillvalue("");
	int partition_count = 1;
        rows = cols = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Parse options
	while ((c = getopt_long(argc, argv, "p:r:f:n:", longopts, NULL)) != -1) {
		switch (c) {
		case 0:
			// getopt_long() set a variable, just keep going
			break;
		case 'p':
			output_srs = optarg;
			break;
		case 'n':
			partition_count = strtol(optarg, NULL, 10);
			break;
		case 'r':
			resampler = optarg;
			break;
		case 'f':
			fillvalue = optarg;
			break;
		default:
			fprintf(stderr, "%s: option '-%c' is invalid: ignored\n", argv[0], optopt);
			break;
		}
	}
	
	if (argc < 2) {
		printf("USAGE %s\n", usage);
		return 0;
	}

	if (driver(argv[argc-2], argv[argc-1], output_srs, fillvalue, partition_count) != 0) {
		MPI_Abort(MPI_COMM_WORLD, 1);
		return 1;
	}
	MPI_Finalize();
	
	return 0;
}
