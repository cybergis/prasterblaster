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

	class RasterChunk 
	{
	public:
		shared_ptr<Projection> projection_;
		Coordinate raster_location;
		Coordinate ul_projected_corner_;
		double pixel_size_; // in meters
                int row_count_;
                int column_count_;
		GDALDataType *pixel_type_;
		int band_count_;
		double geotransform_[6];
		void *pixels_;
	};
	
}

#endif //RASTERCHUNK_HH_
