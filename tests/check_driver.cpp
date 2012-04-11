#include <gtest/gtest.h>

#include <memory>

#include "driver.hh"
#include "projectedraster.hh"
#include "sharedptr.hh"

static string test_dir = "tests/testdata/";
static string output_dir = "tests/testoutput/";

class DriverTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		in = shared_ptr<ProjectedRaster>(new ProjectedRaster(test_dir + "veg_geographic_1deg.tif"));
	}
	shared_ptr<ProjectedRaster> in;

};

TEST(DriverTest, driver_output_creation) {
  
	int ret = driver(test_dir + "veg_geographic_1deg.tif", 
			 output_dir + "veg_mollweide_1deg.tif",
			 "+proj=moll",
			 "max",
			 "",
			 1);
	
	EXPECT_EQ(0, ret);
	
	
	ret = driver(test_dir + "glc_geographic_30sec.tif",
		     output_dir + "glc_mollweide_30sec.tif",
		     "+proj=moll",
		     "max",
		     "",
		     1);
	
	EXPECT_EQ(0, ret);
	
}
