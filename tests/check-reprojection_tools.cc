

#include <gtest/gtest.h>

#include <vector>

#include "src/projectedraster.h"
#include "src/reprojection_tools.h"

using librasterblaster::ProjectedRaster;
using std::vector;

TEST(ReprojectionTest, ProjectedMinbox) {
  ProjectedRaster holdnorm("tests/testdata/holdnorm_geographic_30min.tif");
  ASSERT_TRUE(holdnorm.ready());
}

TEST(ReprojectionTest, RasterMinbox) {
/**
 * In this test we are going to perform some sanity checks on the RasterMinbox function.
 *
 */

  ProjectedRaster holdnorm("tests/testdata/holdnorm_geographic_30min.tif");

  std::string output_srs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs";

  ASSERT_TRUE(holdnorm.ready());

  Coordinate cc(holdnorm.ul_x(),
                holdnorm.ul_y(),
                UNDEF);

  librasterblaster::Area mm = librasterblaster::ProjectedMinbox(cc,
                                                                holdnorm.srs(),
                                                                holdnorm.pixel_size(),
                                                                holdnorm.row_count(),
                                                                holdnorm.column_count(),
                                                                output_srs);

  // Now we will check that mm is well formed
  ASSERT_LT(mm.ul.x, mm.lr.x);
  ASSERT_GT(mm.ul.y, mm.lr.y);
  printf("\t\t%f %f %f %f\n", mm.ul.x, mm.ul.y, mm.lr.x, mm.lr.y);

  // Calculate size of theoretical output raster
  int rows = (mm.ul.y - mm.lr.y) / holdnorm.pixel_size();
  int columns = (mm.lr.x - mm.ul.x) / holdnorm.pixel_size();

  ASSERT_GT(rows, 0);
  ASSERT_GT(columns, 0);
  printf("\t\t%d rows, %d columns\n", rows, columns);


  // Now create a partition on the right edge of the output raster and check its
  // minbox in the input raster.
  librasterblaster::Area partition(0.0, 0.0, columns-1, rows-1);

  Coordinate mmul(mm.ul.x, mm.ul.y, UNDEF);
  shared_ptr<Projection> out_proj(librasterblaster::ProjectionFactory(output_srs));
  librasterblaster::Area minbox = librasterblaster::RasterMinbox(out_proj,
                                                                 mmul,
                                                                 holdnorm.pixel_size(),
                                                                 rows,
                                                                 columns,
                                                                 holdnorm.projection(),
                                                                 cc,
                                                                 holdnorm.pixel_size(),
                                                                 holdnorm.row_count(),
                                                                 holdnorm.column_count(),
                                                                 partition); 

  printf("\t\tDestination raster partition: %f %f %f %f\n", partition.ul.x, partition.ul.y, partition.lr.x, partition.lr.y);
  printf("\t\tSource raster minbox: %f %f %f %f\n", minbox.ul.x, minbox.ul.y, minbox.lr.x, minbox.lr.y);
}

TEST(ReprojectionTest, NA_partitioning) {
  ProjectedRaster nlcd("tests/testdata/nlcd2006_landcover_4-20-11_se5.tif");

  std::string output_srs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs";
  
  ASSERT_TRUE(nlcd.ready());
  
  // In this test we are going to find the minbox in Mollweide and then find the
  // minbox of a partition of the _first_ minbox is well formed.
  
  Coordinate cc(nlcd.ul_x(),
                nlcd.ul_y(),
                UNDEF);

  librasterblaster::Area mm = librasterblaster::ProjectedMinbox(cc,
                                                                nlcd.srs(),
                                                                nlcd.pixel_size(),
                                                                nlcd.row_count(),
                                                                nlcd.column_count(),
                                                                output_srs);
  
  // Now we will check that mm is well formed
  ASSERT_LT(mm.ul.x, mm.lr.x);
  ASSERT_GT(mm.ul.y, mm.lr.y);
  printf("\t\t%f %f %f %f\n", mm.ul.x, mm.ul.y, mm.lr.x, mm.lr.y);

  // Calculate size of theoretical output raster
  int rows = (mm.ul.y - mm.lr.y) / nlcd.pixel_size();
  int columns = (mm.lr.x - mm.ul.x) / nlcd.pixel_size();

  ASSERT_GT(rows, 0);
  ASSERT_GT(columns, 0);
  printf("\t\t%d rows, %d columns\n", rows, columns);

  librasterblaster::Area partition(0, 0, 218920, 0);

  // Create output projection object
  shared_ptr<Projection> out_proj(librasterblaster::ProjectionFactory(output_srs));

  vector<librasterblaster::Area> parts = librasterblaster::RowPartition(0, 1, rows, columns, 5000000);
  printf("\t\t Testing with %lu partitions\n", parts.size());
  Coordinate mmul(mm.ul.x, mm.ul.y, UNDEF);
  librasterblaster::Area minbox;
  for (int i=0 ; i < parts.size(); ++i) {
    partition = parts.at(i);
    //partition.ul.y = rows - partition.ul.y - 1;
    //partition.lr.y = rows - partition.lr.y - 1;
    minbox = librasterblaster::RasterMinbox(out_proj,
                                            mmul,
                                            nlcd.pixel_size(),
                                            rows,
                                            columns,
                                            nlcd.projection(),
                                            cc,
                                            nlcd.pixel_size(),
                                            nlcd.row_count(),
                                            nlcd.column_count(),
                                            partition);
    if ((minbox.ul.x == -1) || (minbox.ul.y == -1) || (minbox.lr.x == -1) || (minbox.lr.y == -1)) {
      printf("\t\t Partition %d", i);
      printf(" partition bad: %f %f %f %f => %f %f %f %f \n", partition.ul.x, partition.ul.y, partition.lr.x, partition.lr.y, minbox.ul.x, minbox.ul.y, minbox.lr.x, minbox.lr.y);


      continue;
    }

    if (minbox.ul.x > minbox.lr.x) {
      printf("\t\tDestination raster partition: %f %f %f %f\n", partition.ul.x, partition.ul.y, partition.lr.x, partition.lr.y);
      printf("\t\tSource raster minbox: %f %f %f %f\n", minbox.ul.x, minbox.ul.y, minbox.lr.x, minbox.lr.y);
    }
    
    EXPECT_LT(minbox.ul.x, minbox.lr.x);
  }
}
