

#include <gtest/gtest.h>

#include <vector>

#include "src/projectedraster.h"
#include "src/reprojection_tools.h"

using librasterblaster::ProjectedRaster;
using std::vector;

TEST(ReprojectionTest, NA_partitioning) {
  ProjectedRaster nlcd("tests/testdata/nlcd2006_landcover_4-20-11_se5.tif");

  std::string output_srs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs";
  
  ASSERT_TRUE(nlcd.ready());
  
  // In this test we are going to find the minbox in Mollweide and then find the
  // minbox of a partition of the _first_ minbox is well formed.
  
  librasterblaster::Area mm = librasterblaster::ProjectedMinbox(Coordinate(nlcd.ul_x(),
                                                                           nlcd.ul_y(),
                                                                           UNDEF),
                                                                nlcd.srs(),
                                                                nlcd.pixel_size(),
                                                                nlcd.row_count(),
                                                                nlcd.column_count(),
                                                                output_srs);
  
  // Now we will check that mm is well formed
  ASSERT_LT(mm.ul.x, mm.lr.x);
  ASSERT_GT(mm.ul.y, mm.lr.y);
  printf("%f %f %f %f\n", mm.ul.x, mm.ul.y, mm.lr.x, mm.lr.y);
  

  // Calculate size of theoretical output raster
  int rows = (mm.ul.y - mm.lr.y) / nlcd.pixel_size();
  int columns = (mm.lr.x - mm.ul.x) / nlcd.pixel_size();

  ASSERT_GT(rows, 0);
  ASSERT_GT(columns, 0);
  printf("%d rows, %d columns\n", rows, columns);

  librasterblaster::Area partition(0, 0, 218920, 0);

  // Create output projection object
  shared_ptr<Projection> out_proj(librasterblaster::ProjectionFactory(output_srs));

  vector<librasterblaster::Area> parts = librasterblaster::RowPartition(0, 1, rows, columns, 2048);
  librasterblaster::Area minbox;
  for (int i=0; i < parts.size(); ++i) {
    partition = parts.at(i);
    //partition.ul.y = rows - partition.ul.y - 1;
    //partition.lr.y = rows - partition.lr.y - 1;
    minbox = librasterblaster::RasterMinbox(nlcd.projection(),
                                            Coordinate(nlcd.ul_x(),
                                                       nlcd.ul_y(),
                                                       UNDEF),
                                            nlcd.pixel_size(),
                                            nlcd.row_count(),
                                            nlcd.column_count(),
                                            out_proj,
                                            Coordinate(mm.ul.x,
                                                       mm.ul.y,
                                                       UNDEF),
                                            nlcd.pixel_size(),
                                            partition);
    if (minbox.ul.x > minbox.lr.x) {
      printf("Destination raster partition: %f %f %f %f\n", partition.ul.x, partition.ul.y, partition.lr.x, partition.lr.y);
      printf("Source raster minbox: %f %f %f %f\n", minbox.ul.x, minbox.ul.y, minbox.lr.x, minbox.lr.y);
    }

    ASSERT_LT(minbox.ul.x, minbox.lr.x);
  }
}
