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
 * The resampler class implements a variety of resampling techniques.
 *
 */


#ifndef RESAMPLER_HH
#define RESAMPLER_HH

#include <gdal.h>

#include "rasterchunk.hh"

namespace Resampler
{	
	template <typename T>
	T Max(Coordinate input_ul,
	      Coordinate input_lr,
	      int input_column_count,
	      T* input_pixels)
	{
		T temp = input_pixels[(int)input_ul.y * input_column_count + (int)input_ul.x];
		T temp2 = 0;
		for (int x = input_ul.x; x <= input_lr.x; ++x) {
			for (int y = input_ul.y; y <= input_lr.y; ++y) {
				temp2 = input_pixels[y * input_column_count + x];
				
				if (temp2 > temp) {
					temp = temp2;
				}
			}
		}
		
		return temp;
	}

	template <typename T>
	T Min(Coordinate input_ul,
	      Coordinate input_lr,
	      int input_column_count,
	      T* input_pixels)
	{
		T temp = input_pixels[(int)input_ul.y * input_column_count + (int)input_ul.x];
		T temp2 = 0;
		for (int x = input_ul.x; x <= input_lr.x; ++x) {
			for (int y = input_ul.y; y <= input_lr.y; ++y) {
				temp2 = input_pixels[y * input_column_count + x];
				
				if (temp2 < temp) {
					temp = temp2;
				}
			}
		}

		return temp;
	}

	template <typename T>
	T Mean(Coordinate input_ul,
	       Coordinate input_lr,
	       int input_column_count,
	       T* input_pixels)
	{
		T temp = 0;
		
		for (int x = input_ul.x; x <= input_lr.x; ++x) {
			for (int y = input_ul.y; y <= input_lr.y; ++y) {
				temp += input_pixels[y * input_column_count + x];
			}
		}

		temp /= (input_lr.x - input_ul.x) * (input_lr.y - input_ul.y);

		return temp;
	}

	template <typename T>
	T Median(Coordinate input_ul,
		 Coordinate input_lr,
		 int input_column_count,
		 T* input_pixels)
	{
		
	}
	
	template <typename T>
	T Mode(Coordinate input_ul,
	       Coordinate input_lr,
	       int input_column_count,
	       T* input_pixels)
	{
		
	}

	template <typename T>
	T Sum(Coordinate input_ul,
		   Coordinate input_lr,
		   int input_column_count,
		   T* input_pixels)
	{
		
	}


	template <typename T>
	T Bilinear(Coordinate input_ul,
		   Coordinate input_lr,
		   int input_column_count,
		   T* input_pixels)
	{
		
	}
} // namespace resampler

#endif // RESAMPLER_HH
