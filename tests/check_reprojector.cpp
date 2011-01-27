/* 
 * Filename: check_reprojector.cpp
 * Author: David Mattli <dmattli@usgs.gov>
 * License: PUBLIC DOMAIN
 */

#include <vector>
#include <algorithm>

#include <gtest/gtest.h>

#include "reprojector.hh"

bool extents_are_continuous(vector<ChunkExtent> extents, long row_count) 
{
	long progress = 0;


	std::sort(extents.begin(), extents.end());

	return true;
}

TEST(reprojector_test, extents_basic_tests) {

	vector<ChunkExtent> ce = Reprojector::getChunkExtents(10000, 100, 10);



	ASSERT_TRUE(extents_are_continuous(ce, 10000));

}
