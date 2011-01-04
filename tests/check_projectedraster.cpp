

#include <gtest/gtest.h>

#include "projectedraster.hh"

TEST(PR_TEST, OpenFile) {
	ProjectedRaster pr("test.tif");
	EXPECT_TRUE(pr.isReady());
}
