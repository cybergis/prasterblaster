/*!
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE 
 *
 * This software is in the public domain, furnished "as is", without
 * technical support, and with no warranty, express or implied, as to
 * its usefulness for any purpose.
 *
 * @section DESCRIPTION
 *
 * 
 *
 */

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
