


#include "rasterchunk.hh"

namespace RasterChunk {

	RasterChunk::~RasterChunk()
	{
		if (this->pixels_ != NULL) {
			free(this->pixels_);
		}

		return;
	}

}
