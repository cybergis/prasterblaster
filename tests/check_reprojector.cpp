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
	vector<ChunkExtent> ce = Reprojector::getChunkExtents(10000, 100, 10);
	
	std::sort(ce.begin(), ce.end());

	ASSERT_EQ(10000 - 1, ce.back().lastIndex());

}

TEST(reprojector_test, extents_continuity) {

	vector<ChunkExtent> ce = Reprojector::getChunkExtents(10000, 100, 10);


	ASSERT_EQ(10000 - 1, ce.back().lastIndex());
	ASSERT_TRUE(extents_are_continuous(ce, 10000));

}
