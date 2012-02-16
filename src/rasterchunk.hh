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
 * The RasterChunk class represents an in-memory, georeferenced section of raster. 
 *
 */


#ifndef RASTERCHUNK_HH_
#define RASTERCHUNK_HH_

#include <memory>

#include <gdal_priv.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/coordinate.h"

using std::shared_ptr;

namespace RasterChunk {

	class RasterChunk 
	{
	public:
		RasterChunk() {};
		~RasterChunk();
		shared_ptr<Projection> projection_;
		Coordinate raster_location_;
		Coordinate ul_projected_corner_;
		double pixel_size_; // in meters
                int row_count_;
                int column_count_;
		GDALDataType pixel_type_;
		int band_count_;
		double geotransform_[6];
		void *pixels_;
	};
	
}

#endif //RASTERCHUNK_HH_
