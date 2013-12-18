/*!
 * Copyright 0000 <Nobody>
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
 * System level tests
 *
 */


#include <vector>

#include <mpi.h>

#include "gtest/gtest.h"
#include "src/configuration.h"
#include "src/quadtree.h"
#include "src/reprojection_tools.h"
#include "src/resampler.h"
#include "src/demos/prasterblaster-pio.h"

#include "src/demos/sptw.h"

using librasterblaster::Configuration;
using librasterblaster::PRB_NOERROR;
using librasterblaster::rastercompare;

namespace {
class MinimalistPrinter : public ::testing::EmptyTestEventListener {
  // Called before a test starts.
  virtual void OnTestStart(const ::testing::TestInfo& test_info) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      printf("*** Test %s.%s starting.\n",
             test_info.test_case_name(), test_info.name());
    }
  }

  // Called after a failed assertion or a SUCCEED() invocation.
  virtual void OnTestPartResult(
      const ::testing::TestPartResult& test_part_result) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      printf("%s in %s:%d\n%s\n",
             test_part_result.failed() ? "*** Failure" : "Success",
             test_part_result.file_name(),
             test_part_result.line_number(),
             test_part_result.summary());
    }
  }

  // Called after a test ends.
  virtual void OnTestEnd(const ::testing::TestInfo& test_info) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      printf("*** Test %s.%s ending.\n",
             test_info.test_case_name(), test_info.name());
    }
  }
};

// The fixture for system testing
class RasterTest : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  RasterTest() {
          // You can do set-up work for each test here.
  }

  virtual ~RasterTest() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Objects Declared here can by used by all tests in fixture
  int rank;
  int process_count;
};

TEST_F(RasterTest, Veg1Deg) {
  Configuration conf;
  conf.partition_size = 1;
  conf.input_filename = "tests/testdata/veg_geographic_1deg.tif";
  conf.output_filename = "tests/testoutput/veg_mollweide_1deg.tif";
  conf.resampler = librasterblaster::MIN;
  conf.output_srs = "+proj=moll +a=6370997 +b=6370997";
  conf.fillvalue = "32";

  int ret = prasterblasterpio(conf);
  ASSERT_EQ(PRB_NOERROR, ret);

  int raster_compare_ret = rastercompare("tests/testdata/veg_mollweide_1deg.tif", 
                                         "tests/testoutput/veg_mollweide_1deg.tif");
  ASSERT_EQ(0, raster_compare_ret);
}

TEST_F(RasterTest, Veg1DegRobin) {
  Configuration conf;
  conf.partition_size = 1;
  conf.input_filename = "tests/testdata/veg_geographic_1deg.tif";
  conf.output_filename = "tests/testoutput/veg_robinson_1deg.tif";
  conf.resampler = librasterblaster::MIN;
  conf.output_srs = "+proj=robin +a=6370997 +b=6370997";
  conf.fillvalue = "32";

  int ret = prasterblasterpio(conf);
}

TEST_F(RasterTest, Veg1DegSinusoidal) {
  Configuration conf;
  conf.partition_size = 1;
  conf.input_filename = "tests/testdata/veg_geographic_1deg.tif";
  conf.output_filename = "tests/testoutput/veg_sinosoidal_1deg.tif";
  conf.resampler = librasterblaster::MIN;
  conf.output_srs = "+proj=sinu +lon_0=0.0 +a=6370997 +b=6370997";
  conf.fillvalue = "32";

  int ret = prasterblasterpio(conf);
}

TEST_F(RasterTest, HoldNorm30Min) {
  Configuration conf;
  conf.partition_size = 1;
  conf.input_filename = "tests/testdata/holdnorm_geographic_30min.tif";
  conf.output_filename = "tests/testoutput/holdnorm_mollweide_30min.tif";
  conf.resampler = librasterblaster::MIN;
  conf.output_srs = "+proj=moll +a=6370997 +b=6370997";
  conf.fillvalue = "32";

  int ret = prasterblasterpio(conf);
  ASSERT_EQ(PRB_NOERROR, ret);

  int raster_compare_ret = rastercompare("tests/testdata/holdnorm_geographic_30min.tif", 
                                         "tests/testoutput/holdnorm_mollweide_30min.tif");
  ASSERT_EQ(0, raster_compare_ret);
}

TEST_F(RasterTest, GLC30sec) {
  Configuration conf;
  conf.partitioner = "tile";
  conf.partition_size = 32;
  conf.input_filename = "tests/testdata/glc_geographic_30sec.tif";
  conf.output_filename = "tests/testoutput/glc_mollweide_30sec.tif";
  conf.resampler = librasterblaster::MIN;
  conf.output_srs = "+proj=moll +a=6370997 +b=6370997";
  conf.fillvalue = "8";

  int ret = prasterblasterpio(conf);
  ASSERT_EQ(PRB_NOERROR, ret);

  int raster_compare_ret = rastercompare("tests/testdata/glc_geographic_30sec.tif", 
                                         "tests/testoutput/glc_geographic_30sec.tif");
  ASSERT_EQ(0, raster_compare_ret);
}

TEST_F(RasterTest, NLCD) {
  Configuration conf;
  conf.partition_size = 64;
  conf.input_filename =
      "tests/testdata/nlcd2006_landcover_4-20-11_se5.tif";
  conf.output_filename =
      "tests/testoutput/nlcd2006_landcover_4-20-11_se5_mollweide.tif";
  conf.resampler = librasterblaster::MIN;
  conf.output_srs =
      "+proj=moll +lat_1=29.5 +lat_2=45.5 "
      "+lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83";
  conf.fillvalue = "32";

  int ret = prasterblasterpio(conf);
}
}  // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
  delete listeners.Release(listeners.default_result_printer());
  listeners.Append(new MinimalistPrinter);
  ::testing::InitGoogleTest(&argc, argv);
  int ret = RUN_ALL_TESTS();
  MPI_Finalize();

  return ret;
}
