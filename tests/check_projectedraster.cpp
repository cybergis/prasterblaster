

#include <gtest/gtest.h>

#include <boost/shared_ptr.hpp>

#include "projectedraster.hh"

using boost::shared_ptr;

static string test_dir = "tests/testdata/";

class ChunkTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		in = shared_ptr<ProjectedRaster>(new ProjectedRaster(test_dir + "veg_geographic_1deg.img"));
	}
	shared_ptr<ProjectedRaster> in;


};


TEST(PR_TEST, OpenFile) {
	ProjectedRaster pr(test_dir + "veg_geographic_1deg.img");
	EXPECT_TRUE(pr.isReady());
}

TEST(PR_TEST, GeoMinbox) {
	ProjectedRaster pr(test_dir + "veg_geographic_1deg.img");
	Area box = pr.getGeographicalMinbox();
	Area box2 = pr.getProjectedMinbox();
	
	printf("Geo Minbox: %f %f %f %f\n", box.ul.x, box.ul.y, box.lr.x, box.lr.y);

}

TEST(PR_TEST, NEW_RASTER) {
	shared_ptr<ProjectedRaster> in(new ProjectedRaster(test_dir + "veg_geographic_1deg.img"));
	shared_ptr<ProjectedRaster> out; 

	ASSERT_TRUE(in->isReady());

	shared_ptr<Projection> in_proj(in->getProjection());
	shared_ptr<Projection> out_proj (Transformer::convertProjection(MOLL));
	out_proj->setUnits(in_proj->units());
	out_proj->setDatum(in_proj->datum());
	out_proj->setParams(in_proj->params());

	bool result = ProjectedRaster::CreateRaster(test_dir + "veg_mollweide_1deg.tif",
						    in,
						    out_proj,
						    in->type,
						    in->pixel_size);

	ASSERT_TRUE(result);

	shared_ptr<ProjectedRaster> newpr(new ProjectedRaster(test_dir + "veg_mollweide_1deg.tif"));

	ASSERT_TRUE(newpr->isReady());

	// Check geographical extent
	Area box = newpr->getGeographicalMinbox();

	printf("New Geo Minbox: %f %f %f %f\n", box.ul.x, box.ul.y, box.lr.x, box.lr.y);

}

TEST_F(ChunkTest, raster_row_continuity) {

	long count = 0;
	
	vector<ChunkExtent> chunks = in->getChunks(10);

	for (int i = 0; i < chunks.size() ; ++i) {
		count += chunks[i].rowCount();
	}

	if (in->getRowCount() != count) {
		printf("Last Chunk: %d first %d last %d size\n", chunks.back().firstIndex(),
		       chunks.back().lastIndex(), chunks.back().rowCount());
	}
	ASSERT_EQ(in->getRowCount(), count);

	
}

TEST_F(ChunkTest, minbox_continuity) {


	vector<ChunkExtent> chunks = in->getChunks(10);
	Area minbox; 

	for (int i = 0; i < chunks.size(); ++i) {
		minbox = chunks[i].getGeographicalMinbox();
		printf("Chunk %d (%f,%f) to (%f,%f)\n", i, 
		       minbox.ul.x, minbox.ul.y,
		       minbox.lr.x, minbox.lr.y);

	}

}
