
#include "configuration.h"
#include "projectedraster.h"
#include "reprojection_tools.h"
#include "sharedptr.h"

#include "sptw.h"

using namespace librasterblaster;
using namespace sptw;

int main(int argc, char *argv[]) {

  // Give MPI_Init first run at the command-line arguments
  MPI_Init(&argc, &argv);
  
  // Initialize Configuration object
  Configuration conf(argc, argv);
  
  if (conf.input_filename == "" || conf.output_filename == "") {
    fprintf(stderr, "Specify an input and output filename\n");
    return 0;
  }
  
  // Open the input raster
  shared_ptr<ProjectedRaster> input_raster(new ProjectedRaster(conf.input_filename));
  if (input_raster->ready() == false) {
    fprintf(stderr, "Error opening input raster!\n");
    return 0;
  }

  // Now we have to create the output raster
  PRB_ERROR err = CreateOutputRaster(input_raster,
                                     conf.output_filename,
                                     input_raster->pixel_size(),
                                     conf.output_srs);

  if (err != PRB_NOERROR) {
    fprintf(stderr, "Error creating raster!\n");
  }
  
  // Open the new output raster
  PTIFF* output_raster = open_raster(conf.output_filename);

  if (output_raster == NULL) {
    fprintf(stderr, "Could not open output raster\n");
    return 0;
  }




  // Clean up
  close_raster(output_raster);
  output_raster = NULL;
  MPI_Finalize();

  return 0;
}
