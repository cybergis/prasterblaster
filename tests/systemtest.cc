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
#include <unistd.h>

#include "gtest/gtest.h"
#include "src/configuration.h"
#include "src/reprojection_tools.h"
#include "src/resampler.h"
#include "src/demos/prasterblaster-pio.h"

#include "src/demos/sptw.h"

#include "tests/rastercompare.h"

using librasterblaster::Configuration;
using librasterblaster::PRB_NOERROR;

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

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

TEST(SystemTest, GLOBALVEG) {
  const int gold_count = 12;
  const std::string golden_rasters[] = { "aea", "cea", "eck4", "eck6", "gall",
                                         "gnom", "laea", "merc", "mill",
                                         "moll", "sinu", "vandg"};
  Configuration conf;

  GDALAllRegister();

  for (unsigned int i = 0; i < gold_count; ++i) {
    const std::string gold_name = STR(__PRB_SRC_DIR__) "/tests/testdata/veg_" + golden_rasters[i] + ".tif";
    const std::string test_name = STR(__PRB_SRC_DIR__) "/tests/testdata/veg_test_" + golden_rasters[i] + ".tif";

    // Open golden raster to extract metadata
    GDALDataset *golden = static_cast<GDALDataset*>(GDALOpen(gold_name.c_str(),
                                                             GA_ReadOnly));

    if (golden == NULL) {
      printf("%s\n", gold_name.c_str());
    }
    ASSERT_TRUE(golden != NULL);

    conf.input_filename = STR(__PRB_SRC_DIR__) "/tests/testdata/veg.tif";
    conf.output_filename = test_name;
    conf.resampler = librasterblaster::MIN;
    conf.output_srs = golden->GetProjectionRef();
    conf.tile_size = 16;
    conf.partition_size = 1;

    GDALClose(golden);

    int ret = prasterblasterpio(conf);
    ASSERT_EQ(PRB_NOERROR, ret);

    int raster_compare_ret = rastercompare(gold_name, test_name);

    ASSERT_EQ(0, raster_compare_ret);
    unlink(test_name.c_str());
  }
}
}  // namespace

int main(int argc, char *argv[]) {
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
