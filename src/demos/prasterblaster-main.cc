
#include <mpi.h>

#include "src/configuration.h"
#include "src/demos/prasterblaster-pio.h"

using namespace librasterblaster;

int main(int argc, char *argv[]) {
  int rank = 0;
  int process_count = 0;
  RasterChunk *in_chunk, *out_chunk;
  
  // Give MPI_Init first run at the command-line arguments
  MPI_Init(&argc, &argv);
  
  // Initialize Configuration object 
  Configuration conf(argc, argv);

  int ret =  prasterblasterpio(conf);
  MPI_Finalize();

  return ret;
}
