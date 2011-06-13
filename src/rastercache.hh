
/**
 * @file
 * @author  David Mattli <dmattli@usgs.gov>
 * @version 1.0
 *
 * @section LICENSE
 *
 * This software was written as part of the official duties of a
 * federal employee and so is in the public domain.
 *
 * @section DESCRIPTION
 *
 * This class provides a memory cached interface to a ProjectedRaster class.
 */

#ifndef RASTERCACHE_HH
#define RASTERCACHE_HH

#include <map>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#include "projectedraster.hh"
#include "reprojector.hh"

using std::string;
using std::vector;


template <class PixelType>
class RasterCache
{
public:
	
	/** 
	 * Constructor that creates a RasterCache by reading the file at the specified path.
	 *
	 * @param filename Filename of raster to open
	 *
	 */
	RasterCache(string filename);
	
	
	/**
	 *
	 * Constructor that creates a RasterCache object from the specified ProjectedRaster object.
	 *
	 * @param raster a shared_ptr to a ProjectedRaster to represent.
	 *
	 */
	RasterCache(shared_ptr<ProjectedRaster> raster);
	~RasterCache();


	/** 
	 * Get the pixel value at the given cartesian point of the raster.
	 *
	 * @param x x coordinate
	 * @param y y coordinate
	 */
	PixelType at(long x, long y);
	PixelType at(int x, int y);

        /** 
	 * Get the pixel value at the given index
	 *
	 * @param i Index of the desired value
	 */
	PixelType at(long i);
	PixelType at(int i);
private:
	bool fetchValue(long i);

	long firstIndex;
	long lastIndex;
	long cacheSize;
	vector<PixelType> pixels;
	shared_ptr<ProjectedRaster> raster;

};
	
template <class PixelType>
RasterCache<PixelType>::RasterCache(string filename)
{
	cacheSize = 5;
	firstIndex = -1;
	lastIndex = -1;

	return;
}

template <class PixelType>
RasterCache<PixelType>::RasterCache(shared_ptr<ProjectedRaster> _raster)
{

	raster = _raster;
	cacheSize = 5;

	return;
}

template <class PixelType>
RasterCache<PixelType>::~RasterCache()
{
	return;
}

template <class PixelType>
PixelType RasterCache<PixelType>::at(long x, long y)
{
	return at(x*y);
}

template <class PixelType>
PixelType RasterCache<PixelType>::at(int x, int y)
{
	return at((long)x, (long)y);
}

template <class PixelType>
PixelType RasterCache<PixelType>::at(long i)
{
	if (i >= firstIndex && i <= lastIndex) {
		try {
			return pixels.at(i-firstIndex);
		} catch (std::out_of_range &e) {
			printf("Tried to Access: %ld, index: %ld vector size: %zd, firstindex: %ld, lastindex %ld\n",
			       i, i-firstIndex, pixels.size(), firstIndex, lastIndex);  
		}
	} else {
		if (!fetchValue(i)) {
			throw std::runtime_error("Error reading input"); 
		}
		return pixels.at(i-firstIndex);
	}
	
	
	return -1;
}


template <class PixelType>
PixelType RasterCache<PixelType>::at(int i)
{
	return at((long)i);
}

template <class PixelType>
bool RasterCache<PixelType>::fetchValue(long i)
{

	long row = i / raster->getColCount();
	long max_index = raster->getRowCount() * raster->getColCount() - 1;
	long rows_to_fetch = cacheSize;
	

	// Determine how many rows to fetch
	if ( i > max_index) {
		stringstream s;
		s << "The index " << i << " is greater than " << max_index;
		firstIndex = -1;
		lastIndex = -1;
		throw std::out_of_range(s.str());
	}

	if ((i + cacheSize * raster->getColCount()) > max_index) {
		rows_to_fetch = raster->getRowCount() - row;
	}

	// Insert try block here and update for other pixel types
	pixels.resize(raster->getColCount() * rows_to_fetch);

	bool result = raster->readRaster(row, rows_to_fetch, &(pixels[0]));

	if (!result) {
		firstIndex = -1;
		lastIndex = -1;
		return false;
	}

	// Update ranges
	firstIndex = i;
	lastIndex = i + rows_to_fetch * raster->getColCount() - 1;
	if (pixels.size() != 0) {
		printf("NOOOOOOOOOOOOOOOOOOOOOOOOO %zd\n\n\n\n\n", pixels.size());
	}

	return true;
}

	

#endif // RASTERCACHE_HH
