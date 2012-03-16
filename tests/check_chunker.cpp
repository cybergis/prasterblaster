



#include <gtest/gtest.h>


#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"

#include "driver.hh"
#include "reprojector.hh"


using std::string;
using boost::shared_ptr;

static string test_dir = "tests/testdata/";

class ChunkerTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		bool result = false;
		in = shared_ptr<ProjectedRaster>(new ProjectedRaster(test_dir + "glc_geographic_30sec.tif"));
		shared_ptr<Projection> proj = shared_ptr<Projection>(Transformer::convertProjection(MOLL));
		proj->setUnits(in->getProjection()->units());
		proj->setDatum(in->getProjection()->datum());
		proj->setParams(in->getProjection()->params());
		ProjectedRaster::CreateRaster(test_dir + "glc_mollweide_30sec-chunk.tif",
					      in, proj, in->getPixelType(), in->getPixelSize());
		out = shared_ptr<ProjectedRaster>(new ProjectedRaster(test_dir + "glc_mollweide_30sec-chunk.tif"));

	}
	
protected:
	shared_ptr<ProjectedRaster> in;
	shared_ptr<ProjectedRaster> out;
};

TEST_F(ChunkerTest, check_fixture) {
	ASSERT_TRUE(out->isReady());
}

TEST_F(ChunkerTest, check_chunk_sizes) {
	Chunker c(in, out);
	vector<ChunkExtent> chunks;
	ASSERT_NO_THROW ({
	    chunks = c.getChunksByCount(10, 1, 0);
		});

	ASSERT_NE(chunks.size(), 0);
	ASSERT_EQ(chunks.back().lastIndex+1, out->getRowCount());
	
}

TEST_F(ChunkerTest, chunk_output_comprehensive) {
	
	Chunker c(in, out);
	RasterCoordTransformer transformer(in, out);
	vector<ChunkExtent> chunks = c.getChunksByCount(5, 10, 2);
	Coordinate temp;


			for (vector<ChunkExtent>::iterator c = chunks.begin();
			     c != chunks.end(); c++) {
				for (int x = 0; x < out->getColCount(); ++x) {
					for (int y = c->firstIndex; y <= c->lastIndex; ++y) {
						temp.x = x;
						temp.y = y;
						try {
						    transformer.Transform(temp);
						} catch (std::exception) {
						  // ignore this point
						  continue;
						}
					}
				}
			}

		
		
}	
