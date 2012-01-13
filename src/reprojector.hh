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
using RasterChunk::RasterChunk;

class RasterCoordTransformer
{
public:
	RasterCoordTransformer(shared_ptr<ProjectedRaster> source, shared_ptr<ProjectedRaster> dest);
	RasterCoordTransformer(shared_ptr<ProjectedRaster> source,
			       shared_ptr<Projection> destination_projection,
			        Coordinate destination_ul,
			       double destination_pixel_size);
	RasterCoordTransformer(shared_ptr<Projection> source_projection,
			       Coordinate source_ul,
			       double source_pixel_size,
			       shared_ptr<Projection> destination_projection,
			       Coordinate destination_ul,
			       double destination_pixel_size);
	~RasterCoordTransformer();
	Area Transform(Coordinate source);
private:
	shared_ptr<Projection> src_proj, dest_proj;
	Coordinate source_ul_;
	double source_pixel_size_;
	Coordinate destination_ul_;
	double destination_pixel_size_;
};

vector<Area> ParititionByCount(shared_ptr<ProjectedRaster> destination,
		       int partition_count);

Area FindOutputArea(shared_ptr<ProjectedRaster> input,
		    shared_ptr<Projection> output_projection,
		    double output_pixel_size);

Area MapDestinationAreatoSource(shared_ptr<ProjectedRaster> source,
				shared_ptr<Projection> destination_projection,
				Area destination_area,
				double destination_pixel_size);

ProjectedRaster* GetOutputRaster(ProjectedRaster* input,
				 Projection *out_proj,
				 double out_pixel_size);

#endif // REPROJECTOR_HH
