/*!
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE
 * This work was produced as a part of the official duties of a
 * federal employee and is in the public domain.

 * @section DESCRIPTION
 *
 *
 */


#ifndef RESAMPLER_HH
#define RESAMPLER_HH

#include <gdal.h>

#include "rasterchunk.hh"

using RasterChunk::RasterChunk;

namespace resampler 
{
/*
	void NearestNeighbor(RasterChunk *source,
			     RasterChunk *dest,
			     int dest_pixel_index,
			     int source_box_indices[4],
			     int band_count)
	{
		dest->pixels[dest_pixel_index] = source->pixels[source_box_indices[0]];
		
		return;
	}
	
	
	
	void bilinear(void *inraster, double in_x,
		      double in_y, unsigned long in_cols,
		      void *outraster,
		      unsigned long out_x, unsigned long out_y,
		      unsigned long out_cols);
	
	
*/	
	
	
} // namespace resampler

#endif // RESAMPLER_HH
