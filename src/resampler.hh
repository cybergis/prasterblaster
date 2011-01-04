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

namespace resampler 
{
	enum RESAMPLERS {
		NO_RESAMPLER = 0,
		NEAREST_NEIGHBOR,
		BILINEAR
	};

	typedef void(*resampler_func)(void *src_raster,
				      long center_index,
				      long pixel_width,
				      long pixel_height,
				      long *index_array,
				      void *dest_pixel);

	resampler_func SelectResampler333(RESAMPLERS method, GDALDataType type);

		

	template <typename T>
	void nearest_neighbor(void *src_raster,
			      long center_index,
			      long pixel_width,
			      long pixel_height,
			      long *index_array,
			      void *dest_pixel)
	{
		long temp = pixel_width = pixel_height = *index_array; // remove this line
		temp++;

		*(T*)dest_pixel = ((T*)src_raster)[center_index];
		return;
	}

	
	
	template <typename T>
	void bilinear(void *inraster, double in_x,
		      double in_y, unsigned long in_cols,
		      void *outraster,
		      unsigned long out_x, unsigned long out_y,
		      unsigned long out_cols);




	
} // namespace resampler

#endif // RESAMPLER_HH
