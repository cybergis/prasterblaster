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

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mollweide.h"

#include "driver.hh"
#include "reprojector.hh"

using std::string;
using std::vector;

static string test_dir = "tests/testdata/";

class ReprojTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		in = shared_ptr<ProjectedRaster>(new ProjectedRaster(test_dir + "/glc_geographic_30sec.img"));
	}
	shared_ptr<ProjectedRaster> in;
};

TEST_F(ReprojTest, paritition_test)
{
	ASSERT_TRUE(in->isReady());

	vector<Area> parts = PartitionByCount(in, 10);

	ASSERT_TRUE(parts.size() == 10);

	EXPECT_EQ(parts.at(0).ul.x, 0);

	EXPECT_EQ(in->getColCount() - 1, parts.back().lr.x);	

}

TEST_F(ReprojTest, parition_comprehensiveness)
{
	vector<Area> parts = PartitionByCount(in, 10);

	ASSERT_EQ(0, parts.at(0).ul.x);
	ASSERT_EQ(0, parts.at(0).ul.y);

	for (int i = 0; i < parts.size() - 1; ++i) {
		ASSERT_EQ(parts.at(i).lr.y, parts.at(i+1).ul.y - 1);

		ASSERT_EQ(0, parts.at(i).ul.x);
		ASSERT_EQ(in->getColCount() - 1, parts.at(i).lr.x);

	}
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	::testing::InitGoogleTest(&argc, argv);

	if (argc > 1) {
		test_dir = argv[1];
	}



	return RUN_ALL_TESTS();
	MPI_Finalize();
	return 0;
}

