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

#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include <gdal.h>

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

vector<Area> PartitionByCount(shared_ptr<ProjectedRaster> destination,
		       int partition_count);

Area FindOutputArea(shared_ptr<ProjectedRaster> input,
		    shared_ptr<Projection> output_projection,
		    double output_pixel_size);

Area MapDestinationAreatoSource(shared_ptr<ProjectedRaster> source,
				shared_ptr<Projection> destination_projection,
				Coordinate destination_ul_corner,
				double destination_pixel_size,
				Area destination_raster_area);

bool ReprojectChunk(RasterChunk::RasterChunk source, RasterChunk::RasterChunk destination);

template <class pixelType>
bool ReprojectChunkType(RasterChunk::RasterChunk source, RasterChunk::RasterChunk destination)
{

	shared_ptr<Projection> outproj, inproj;
        Coordinate temp1, temp2;
	std::vector<char> inraster, outraster;

        outproj = destination.projection_;
        inproj = source.projection_;
	
	Area pixelArea;
	int count = 0;
	int total = 0;

	RasterCoordTransformer rt(destination.projection_, 
				  destination.ul_projected_corner_,
				  destination.pixel_size_,
				  source.projection_,
				  source.ul_projected_corner_,
				  source.pixel_size_);        


	for (int chunk_y = destination.row_count_-10; chunk_y < destination.row_count_; ++chunk_y) {
		for (int chunk_x = 0; chunk_x < destination.column_count_; ++chunk_x) {
			
		}

	}

	for (int chunk_y = 0; chunk_y < destination.row_count_; ++chunk_y)  {
		for (int chunk_x = 0; chunk_x < destination.column_count_; ++chunk_x) {
			temp1.x = chunk_x; 
			temp1.y = chunk_y;
			try {
				pixelArea = rt.Transform(temp1);   /// Now makes sense.
			} catch (std::runtime_error) {
				continue;
			}

			if (pixelArea.ul.x == -1.0) {
			  continue;
			}

			temp1 = pixelArea.ul;
			temp2 = pixelArea.lr;

			long ul_x = (long)temp1.x;
			long ul_y = (long)temp1.y;
			long lr_x = (long)temp2.x;
			long lr_y = (long)temp2.y;

			if (ul_x < 0) {
				ul_x = 0;
			}
			
			if (ul_y < 0) {
				ul_y = 0;
			}

			if (lr_x > (source.column_count_ - 1)) {
				lr_x = source.column_count_ - 1;
			}

			if (lr_y > (source.row_count_ - 1)) {
				lr_y = source.row_count_ - 1;
			}
			

			// TODO: Check that ul/lr actually enclose an area
			//       Otherwise, use nearest-neighbor

			// Perform resampling...
			// Write pixel to destination
			reinterpret_cast<pixelType*>(destination.pixels_)[chunk_x + chunk_y * destination.column_count_] = 
				reinterpret_cast<pixelType*>(source.pixels_)[ul_x + ul_y * source.column_count_];

		}
	}

	return true;
}

ProjectedRaster* GetOutputRaster(ProjectedRaster* input,
				 Projection *out_proj,
				 double out_pixel_size);

#endif // REPROJECTOR_HH
