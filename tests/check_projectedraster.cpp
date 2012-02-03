

#include <gtest/gtest.h>

#include <memory>

#include "projectedraster.hh"

using std::shared_ptr;

static string test_dir = "tests/testdata/";
static string output_dir = "tests/testoutput/";

class ChunkTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		in = shared_ptr<ProjectedRaster>(new ProjectedRaster(test_dir + "veg_geographic_1deg.img"));
		
		shared_ptr<Projection> outproj(Transformer::convertProjection(MOLL));
		outproj->setUnits(METER);
		outproj->setDatum(in->getDatum());
		ProjectedRaster::CreateRaster(output_dir + "veg_mollweide_1deg_test.tif",
					      in,
					      outproj,
					      in->getPixelType(),
					      in->getPixelSize());
		out = shared_ptr<ProjectedRaster>(new ProjectedRaster(output_dir + "veg_mollweide_1deg_test.tif"));
	
								   
								      
								      
	}
	shared_ptr<ProjectedRaster> in;
	shared_ptr<ProjectedRaster> out;

};

TEST_F(ChunkTest, raster_row_continuity) {

  /*	long count = 0;
	
	vector<ChunkExtent> chunks = in->getChunks(10);

	for (int i = 0; i < chunks.size() ; ++i) {
		count += chunks[i].rowCount();
	}

	if (in->getRowCount() != count) {
		printf("Last Chunk: %d first %d last %d size\n", chunks.back().firstIndex(),
		       chunks.back().lastIndex(), chunks.back().rowCount());
	}
	ASSERT_EQ(in->getRowCount(), count);
  */
	
}

TEST_F(ChunkTest, raster_writing) {
	

}

TEST_F(ChunkTest, minbox_continuity) {

	RasterChunk::RasterChunk *chunk = NULL;
	Area area;

	area.ul.x = 0;
	area.ul.y = 0;
	area.lr.x = in->getColCount() - 1;
	area.lr.y = in->getRowCount() - 1;

	chunk = in->createRasterChunk(area);

	ASSERT_NE(chunk, (void*)NULL);

	EXPECT_EQ(in->getPixelType(), chunk->pixel_type_);

	delete chunk;

}
