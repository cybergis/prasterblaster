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
 *
 */


#ifndef RASTERCHUNK_HH_
#define RASTERCHUNK_HH_


#include <gdal_priv.h>


#include "projectedraster.hh"

namespace RasterChunk {

        template<type t> class Coordinate
        {
        public:
                t x;
                t y;
        };

	class RasterChunk 
	{
	public:
		shared_ptr<Projection> projection;
                Coordinate<int> upper_left_location;
                int row_count;
                int column_count;
		GDALDataType *pixel_type;
		int band_count;
		double geotransform[6];
		void *pixels;
	};
	
}

#endif //RASTERCHUNK_HH_
