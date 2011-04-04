



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
		in = shared_ptr<ProjectedRaster>(new ProjectedRaster(test_dir + "glc_geographic_30sec.img"));
		shared_ptr<Projection> proj = shared_ptr<Projection>(Transformer::convertProjection(MOLL));
		proj->setParams(in->getProjection()->params());
		out = shared_ptr<ProjectedRaster>(new ProjectedRaster(test_dir + "glc_mollweide_30sec-chunk.img"));
	}
	
protected:
	shared_ptr<ProjectedRaster> in;
	shared_ptr<ProjectedRaster> out;
};


TEST_F(ChunkerTest, chunk_reasonableness) {
	
	Chunker c(in, out);
	
	vector<ChunkExtent> chunks = c.getChunksByCount(10, 20);
	
}	
