/* 
 * Filename: check_reprojector.cpp
 * Author: David Mattli <dmattli@usgs.gov>
 * License: PUBLIC DOMAIN
 */


#include <algorithm>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <gtest/gtest.h>
#include <mpi.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mollweide.h"

#include "driver.hh"
#include "reprojector.hh"

using std::string;
using boost::shared_ptr;

static string test_dir;

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


TEST(reprojector_test, minbox_test) {
	
	string input_raster = test_dir + "/testdata/veg_geographic_1deg.img";
	string output_raster = test_dir + "/testdata/minbox_veg_geographic_1deg.tif";
	shared_ptr<ProjectedRaster> in(new ProjectedRaster(input_raster));
	shared_ptr<Projection> in_proj, out_proj(Transformer::convertProjection(MOLL));

	ASSERT_TRUE(in->isReady());
	in_proj = shared_ptr<Projection>(in->getProjection());
	
	out_proj->setUnits(in_proj->units());
	out_proj->setDatum(in_proj->datum());
	out_proj->setParams(in_proj->params());

	bool result = ProjectedRaster::CreateRaster(output_raster,
						    in,
						    out_proj,
						    in->type,
						    in->pixel_size);

	ASSERT_TRUE(result);

	shared_ptr<ProjectedRaster> out(new ProjectedRaster(output_raster));
	
	ASSERT_TRUE(out->isReady());

	

	

	
	

}

TEST(reprojector_test, mollweide_projection_test) {
	
	string input_raster = test_dir + "/testdata/veg_geographic_1deg.img";
	string output_raster = test_dir + "/testdata/veg_mollweide_1deg.tif"; 


	printf("Reprojecting: %s\n", input_raster.c_str());
	int ret = driver(input_raster, output_raster, "mollweide");

	ASSERT_EQ(0, ret);

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

