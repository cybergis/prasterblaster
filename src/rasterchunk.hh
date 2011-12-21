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

	class Interval
	{
	public:
		Interval() {};
		Interval(int _first_index, int _last_index) : first_index(_first_index),
                                                                          last_index(_last_index) {}
		int first_index;
 		int last_index;
	};
	
	class IntervalPair
	{
	public:
		IntervalPair() {};
		IntervalPair(Interval _source, Interval _destination) : source(_source),
									destination(_destination) {}
									       
		Interval getSourceInterval() { return source; }
		Interval getDestinationInterval() { return destination; }
		void setSourceInterval(Interval _source) { source = _source; }
		void setDestinationInterval(Interval _destination) { destination = _destination; }
	     
	private:
		Interval source;
		Interval destination;
	};

	class RasterChunk 
	{
	public:
		shared_ptr<Projection> projection;
		Interval extent;
		GDALDataType *pixel_type;
		int band_count;
		double geotransform[6];
		void *pixels;
	};
	
}

#endif //RASTERCHUNK_HH_
