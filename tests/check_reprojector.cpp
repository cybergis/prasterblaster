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

static string test_dir = "tests/";

class ReprojTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		in = shared_ptr<ProjectedRaster>(new ProjectedRaster(test_dir + "glc_geographic_30sec.img"));
	}
	shared_ptr<ProjectedRaster> in;
};


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

	ASSERT_FALSE(out_proj->errorOccured());

	bool result = ProjectedRaster::CreateRaster(output_raster,
						    in,
						    out_proj,
						    in->type,
						    in->pixel_size);

	ASSERT_FALSE(out_proj->errorOccured());

	ASSERT_TRUE(result);

	shared_ptr<ProjectedRaster> out(new ProjectedRaster(output_raster));
	
	ASSERT_TRUE(out->isReady());

}

TEST(reprojector_test, mollweide_projection_test) {
	
	string input_raster = test_dir + "/testdata/veg_geographic_1deg.img";
	string output_raster = test_dir + "/testdata/veg_mollweide_1deg.tif"; 
	int ret = 0;

	  printf("Reprojecting: %s\n", input_raster.c_str());

	ASSERT_NO_THROW ({
	    try {
	    ret = driver(input_raster, output_raster, "mollweide");
	    } catch(std::exception &e) {
	      printf("Error! ");
	      std::cout << e.what();
	      throw e;
	    }
	  });
	ASSERT_EQ(0, ret);

}

TEST(reprojector_test, chunk_mollweide) {
	
	string input_raster = test_dir + "/testdata/glc_geographic_30sec.img";
	string output_raster = test_dir + "/testdata/glc_mollweide_30sec.tif"; 
	
	
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

