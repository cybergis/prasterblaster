/* 
 * Filename: check_reprojector.cpp
 * Author: David Mattli <dmattli@usgs.gov>
 * License: PUBLIC DOMAIN
 */


#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <mpi.h>

#include "driver.hh"
#include "reprojector.hh"

using std::string;

string test_dir;

bool extents_are_continuous(vector<ChunkExtent> extents, long row_count) 
{
	long progress = 0;


	std::sort(extents.begin(), extents.end());

	for (int i = 0; i < extents.size()-1; ++i) {
		if (extents[i].lastIndex() +1 != extents[i+1].firstIndex()) {
			fprintf(stderr, "Discontinuity at index: %ld, values: %ld %ld\n",
				i, extents[i].lastIndex()+1, extents[i+1].firstIndex());
			return false;
		}
	}

	return true;
}

TEST(reprojector_test, extents_correct_max_value) {
	
	for (int i = 1000; i < 20000; i+=1000) {
		for (int j=357; j < 500; j+=100) {
			for (int k=j; k < j*2; k+=100) {
				vector<ChunkExtent> ce = Reprojector::getChunkExtents(i, k, j);
				std::sort(ce.begin(), ce.end());
				ASSERT_EQ(i - 1, ce.back().lastIndex());
				
			}
		}
		
	}
}

TEST(reprojector_test, extents_continuity) {


	vector<ChunkExtent> ce = Reprojector::getChunkExtents(10000, 100, 10);

	for (int i = 1; i <= 10000; i++) {
		ce = Reprojector::getChunkExtents(i, 100, 10);
		ASSERT_EQ(i - 1, ce.back().lastIndex());
		

	}

	ASSERT_EQ(10000 - 1, ce.back().lastIndex());
	ASSERT_TRUE(extents_are_continuous(ce, 10000));

}


TEST(reprojector_test, basic_reprojection) {
	
	string input_raster = test_dir + "/testdata/veg_geographic_1deg.img";
	string output_raster = test_dir + "/testdata/veg_mollweide_1deg.tif"; 
	string output_raster_desc = test_dir + "/testdata/veg_mollweide_1deg.xml";

	printf("Reprojecting: %s\n", input_raster.c_str());
	int ret = driver(input_raster, output_raster, output_raster_desc);

	ASSERT_EQ(ret, 0);

}


int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	test_dir = "tests";

	::testing::InitGoogleTest(&argc, argv);

	if (argc > 1) {
		test_dir = argv[1];
	}



	return RUN_ALL_TESTS();
	MPI_Finalize();
	return 0;
}

