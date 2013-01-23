
#include <gtest/gtest.h>

#include <vector>

#include "src/projectedraster.h"
#include "src/reprojection_tools.h"
#include "src/rastercoordtransformer.h"
#include "src/gctp_cpp/transformer.h"

using librasterblaster::ProjectedRaster;
using librasterblaster::RasterCoordTransformer;
using librasterblaster::Area;
using std::vector;

TEST(RasterCoordTransformer, EdgeTransformations) {
  shared_ptr<ProjectedRaster> veg(new ProjectedRaster("tests/testdata/veg_geographic_1deg.tif"));
  
  ASSERT_TRUE(veg->ready());
  
  shared_ptr<Projection> oldproj(veg->projection());
  shared_ptr<Projection> newproj(Transformer::convertProjection(MOLL));

  newproj->setUnits(METER);
  newproj->setDatum(WGS_84);
  
  Coordinate new_ul;
  oldproj->inverse(veg->ul_x(), veg->ul_y(), &new_ul.x, &new_ul.y);
  printf("Geographical UL: %f %f\n", new_ul.x, new_ul.y);
  newproj->forward(new_ul.x, new_ul.y, &new_ul.x, &new_ul.y);
  printf("New projected UL: %f %f\n", new_ul.x, new_ul.y);
  
  RasterCoordTransformer rct(veg,
                             newproj,
                             new_ul,
                             veg->pixel_size());
  
  Coordinate coord;
  coord.x = veg->column_count() - 1;
  coord.y = veg->row_count() - 1;
  printf("%f %f maps to:\n", coord.x, coord.y);
  Area a = rct.Transform(coord);
  

  printf("Upper-left: %f %f ... Lower-right: %f %f \n", a.ul.x, a.ul.y, a.lr.x, a.lr.y);

  coord.x = 0.0;
  coord.y = 0.0;
  printf("%f %f maps to:\n", coord.x, coord.y);
  a = rct.Transform(coord);
  printf("Upper-left: %f %f ... Lower-right: %f %f \n", a.ul.x, a.ul.y, a.lr.x, a.lr.y);

  return;
}
