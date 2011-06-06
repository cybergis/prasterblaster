

/* 
 * Filename: check_rastercache.cpp
 * Author: David Mattli <dmattli@usgs.gov>
 * License: PUBLIC DOMAIN
 */

#include <vector>

#include <boost/shared_ptr.hpp>
#include <gtest/gtest.h>

#include "rastercache.hh"
#include "driver.hh"


using std::vector;
using boost::shared_ptr;


static string test_dir = "tests/testdata/";

class CacheTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		small_raster = shared_ptr<ProjectedRaster>(
			new ProjectedRaster(test_dir + "veg_geographic_1deg.img"));
	}

	shared_ptr<ProjectedRaster> small_raster;

};


TEST_F(CacheTest, fixture) {
	// Check that fixture is in good state
	if (!small_raster->isReady()) {
		fprintf(stderr, "Error opening file \"%s\"in fixture\n", small_raster->filename.c_str());
	}

	ASSERT_TRUE(small_raster->isReady());



}

TEST_F(CacheTest, allocation) {
	RasterCache<unsigned char> cache(small_raster);

	long rows = small_raster->getRowCount();
	long cols = small_raster->getColCount();
	unsigned char t = 0;

	for (int i = 0; i < rows*cols; ++i) {
		ASSERT_NO_THROW({
			    try {
			    t = cache.at(i);
			    } catch (std::exception &e) {
				    fprintf(stderr, "Exception caught! %s %d\n", e.what(), i);
				    throw e;
			    }
		    });
	}


}


TEST_F(CacheTest, pixel_values) {
	RasterCache<unsigned char> cache(small_raster);
	vector<unsigned char> temp_row(small_raster->getColCount());
	unsigned int a(0), b(0);

	ASSERT_NO_THROW ({
			for (int i = 0; 
			     i < small_raster->getColCount() * small_raster->getRowCount();
			     ++i) {
				a = cache.at(i);
				long row = i / small_raster->getColCount();
				long col = i % small_raster->getColCount();
				small_raster->readRaster(row, 1, &(temp_row[0]));
				b = temp_row[col];
				
				ASSERT_EQ(a, b);
			}
		});

}

TEST_F(CacheTest, random_values) {



}
