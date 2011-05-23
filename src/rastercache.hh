

#ifndef RASTERCACHE_HH
#define RASTERCACHE_HH

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "projectedraster.hh"
#include "reprojector.hh"

using std::string;
using std::vector;


template <class PixelType>
class RasterCache
{
public:
	RasterCache(string filename);
	RasterCache(shared_ptr<ProjectedRaster> raster);
	~RasterCache();

	PixelType at(long x, long y);
	PixelType at(int x, int y);
	PixelType at(long i);
	PixelType at(int i);

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
	PixelType t;

	if (i >= firstIndex && i <= lastIndex) {
		return pixels.at(i-firstIndex);
	} else {
		fetchValue(i);
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
	long col = i % raster->getColCount();
	long max_index = raster->getRowCount() * raster->getColCount() - 1;
	long rows_to_fetch = cacheSize;
	

	// Determine how many rows to fetch
	if ( i > max_index) {
		throw std::out_of_range("");
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

	return true;
}

	



#endif // RASTERCACHE_HH
