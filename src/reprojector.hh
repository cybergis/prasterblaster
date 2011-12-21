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

/*! \mainpage Reprojector class for pRasterBlaster
 *
 * \section desc_sec Description 
 * The Reprojector class takes two
 * rasters and reprojects and resamples the source to fill the destination.
 *
 *
 */

#ifndef REPROJECTOR_HH
#define REPROJECTOR_HH

#include <memory>
#include <vector>

#include "gctp_cpp/projection.h"

#include "projectedraster.hh"
#include "rasterchunk.hh"
#include "resampler.hh"

using std::shared_ptr;
using RasterChunk::IntervalPair;
using RasterChunk::Interval;
using RasterChunk::RasterChunk;

class RasterCoordTransformer
{
public:
	RasterCoordTransformer(shared_ptr<ProjectedRaster> source, shared_ptr<ProjectedRaster> dest);
	~RasterCoordTransformer();
	Area Transform(Coordinate source);
private:
	shared_ptr<ProjectedRaster> src, dest;
	shared_ptr<Projection> src_proj, dest_proj;
};


vector<Interval> ParititionByCount(shared_ptr<ProjectedRaster> destination,
				   int partition_count);

Interval FindSourceInterval(shared_ptr<ProjectedRaster> source,
			    shared_ptr<ProjectedRaster> destination,
			    Interval dest_area);
/*
bool ReprojectChunk(RasterChunk *source, RasterChunk *destination);
*/
Area FindRasterArea(shared_ptr<ProjectedRaster> source_raster,
		    shared_ptr<ProjectedRaster> dest_raster,
		    int first_row,
		    int last_row);


Area FindGeographicalExtent(shared_ptr<Projection> projection, 
			    Coordinate ul_point,
			    int row_count,
			    int col_count,
			    double pixel_size);

Area FindProjectedExtent(shared_ptr<Projection> projection,
			 Area geographical_area,
			 double pixel_size);


/*!
 *
 */
Area FindMinBox(double in_ul_x, double in_ul_y,
		 double in_pix_size,
		 int rows, int cols, 
		 Projection *inproj,
		 Projection *outproj,
		 double out_pixsize);


ProjectedRaster* GetOutputRaster(ProjectedRaster* input,
				 Projection *out_proj,
				 double out_pixel_size);

#endif // REPROJECTOR_HH
