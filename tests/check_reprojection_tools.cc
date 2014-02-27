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
 * Tests over the reprojection tools
 *
 */

#include <gtest/gtest.h>

#include <vector>

#include "src/utils.h"
#include "src/reprojection_tools.h"

using librasterblaster::Area;
using librasterblaster::PartitionTile;
using std::vector;

TEST(PartitionTile, SmallRasterManyProcesses) {
  const int process_count = 1000;
  const int64_t row_count = 180;
  const int64_t column_count = 360;
  const int64_t tile_size = 1024;
  const int64_t partition_size = 1;
  std::vector<Area> p = PartitionTile(0,
                                      process_count,
                                      row_count,
                                      column_count,
                                      tile_size,
                                      tile_size,
                                      partition_size);

  ASSERT_EQ(1, p.size());

  p = PartitionTile(5,
                    process_count,
                    row_count,
                    column_count,
                    tile_size,
                    tile_size,
                    partition_size);
}

TEST(PartitionTile, BigRasterManyProcesses) {
  const int process_count = 1000;
  const int64_t row_count = 21600;
  const int64_t column_count = 43200;
  const int64_t tile_size = 1024;
  const int64_t partition_size = 1;
  std::vector<Area> p = PartitionTile(0,
                                      process_count,
                                      row_count,
                                      column_count,
                                      tile_size,
                                      tile_size,
                                      partition_size);
  ASSERT_EQ(1, p.size());
}

TEST(PartitionTile, HugeRasterManyProcesses) {
  const int process_count = 1000;
  const int64_t row_count = 104424;
  const int64_t column_count = 161190;
  const int64_t tile_size = 1024;
  const int64_t partition_size = 1;
  std::vector<Area> p = PartitionTile(0,
                                      process_count,
                                      row_count,
                                      column_count,
                                      tile_size,
                                      tile_size,
                                      partition_size);
  ASSERT_NE(1, p.size());

}

