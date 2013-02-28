
#include <mpi.h>

#include <gtest/gtest.h>
#include "src/configuration.h"
#include "src/projectedraster.h"
#include "src/reprojection_tools.h"
#include "src/demos/prasterblaster-pio.h"

using librasterblaster::Configuration;
using std::string;

const char *input_rasters_[] = { "glc_geographic_30sec.tif",
                                 "holdnorm_geographic_30min.tif",
                                 "nlcd2006_landcover_4-20-11_se5.tif",
                                 "veg_geographic_1deg.tif" };


class SystemTest : public ::testing::Test {
 public:
  const char *output_projection_ = "+proj=moll";
  int rank_;
  int process_count_;
  
 protected:
  virtual void SetUp() {
    MPI_Init(0, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_count_);

  }
  
  virtual void TearDown() {
    MPI_Finalize();
  }
  
};

TEST_F(SystemTest, glc_geographic_30sec_tif) {
  Configuration conf;

  conf.input_filename = string("tests/testdata/") + input_rasters_[0];
  conf.output_filename = "tests/testoutput/glc_mollweide_30sec.tif";
  conf.output_srs = output_projection_;

  prasterblaster_main(conf, rank_, process_count_);

}
