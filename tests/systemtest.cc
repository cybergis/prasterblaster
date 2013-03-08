

#include <vector>

#include "src/configuration.h"
#include "src/projectedraster.h"
#include "src/quadtree.h"
#include "src/reprojection_tools.h"
#include "src/resampler.h"
#include "src/sharedptr.h"
#include "src/demos/prasterblaster-pio.h"

#include "src/demos/sptw.h"

using librasterblaster::Configuration;

int main(int argc, char *argv[]) {
  const int raster_count = 4;
  const char *input_files[] = { 
    "tests/testdata/veg_geographic_1deg.tif",
    "tests/testdata/holdnorm_geographic_30min.tif",
    "tests/testdata/glc_geographic_30sec.tif",
    "tests/testdata/nlcd2006_landcover_4-20-11_se5.tif" 
  };

  const char *output_files[] = {
    "tests/testoutput/veg_geographic_1deg.tif",
    "tests/testoutput/holdnorm_geographic_30min.tif",
    "tests/testoutput/glc_mollweide_30sec.tif",
    "tests/testoutput/nlcd2006_landcover_4-20-11_se5_mollweide.tif"
  };

  const char *output_srs[] = {
    "+proj=moll +a=6370997 +b=6370997",
    "+proj=moll +a=6370997 +b=6370997",
    "+proj=moll +a=6370997 +b=6370997",
    "+proj=moll +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83"
  };

  int rank = 0;
  int process_count = 0;
  double runtimes[raster_count] = {0.0};

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);

  Configuration conf(argc, argv);
  if (conf.partition_size == -1) {
    conf.partition_size = 21600;
  }
  
  conf.resampler = librasterblaster::MIN;

  int job_size = 2;
  if ((argc > 1) && strcmp(argv[argc-1], "full") == 0) {
    job_size = raster_count;
  }

  int ret = 0;
  double start_time, stop_time;
  for (int i=0; i<job_size; ++i) {
    conf.input_filename = input_files[i];
    conf.output_filename = output_files[i];
    conf.output_srs = output_srs[i];

    start_time = MPI_Wtime();
    ret = prasterblaster_main(conf, rank, process_count);
    stop_time = MPI_Wtime();
    runtimes[i] = stop_time - start_time;

    if (ret != 0) {
      fprintf(stderr, "Error reprojecting: %s\n", conf.input_filename.c_str());
      MPI_Finalize();
      return 1;
    }

  }
  // Now print runtimes for each raster
  if (rank == 0) {
          printf("\n\nRUNTIMES:\n\n");
          for (int i=0; i<job_size; ++i) {
                  printf("FILE:\t%s", output_files[i]);
                  printf("\t%f seconds\n", runtimes[i]);
          }
  }

  MPI_Finalize();

  return 0;
}
