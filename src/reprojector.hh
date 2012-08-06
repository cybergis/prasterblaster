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
 * This file implements the main reprojection functionality.
 *
 */

/*! \mainpage Reprojector class for pRasterBlaster
 *
 *
 */

#ifndef REPROJECTOR_HH
#define REPROJECTOR_HH

#include <limits>
#include <stdexcept>
#include <vector>

#include <gdal.h>

#include "gctp_cpp/projection.h"

#include "projectedraster.hh"
#include "rasterchunk.hh"
#include "resampler.hh"
#include "sharedptr.hh"

//! An enum 
/*!
 * Enumerates the types of resamplers available.
 */
enum RESAMPLER {
	NEAREST,
	MIN,
	MAX,
	MEAN,
	MEDIAN,
	MODE,
	SUM, 
	BILINEAR
};


enum PRB_ERROR {
	NO_ERROR,
	FILE_ERROR,
	PROJ_ERROR
};

//! Raster Coordinate transformation class */
/*!
  This class implements the transformation of raster coordinates between two raster spaces with different projections and scales.
 */
class RasterCoordTransformer
{
public:
	//! A constructor
	/*!  

	  This constructor takes two initialized ProjectedRaster
	  shared_ptrs and constructs a ready RasterCoordTransformer.
	 */
	RasterCoordTransformer(shared_ptr<ProjectedRaster> source, shared_ptr<ProjectedRaster> dest);

        //! A constructor
	/*!

	  This constructor takes a shared_ptr to a source raster and
	  three parameters that describe a destination raster: a
	  projection, upper-left coordinate, and a pixel size in
	  meters.

	 */
	RasterCoordTransformer(shared_ptr<ProjectedRaster> source,
			       shared_ptr<Projection> destination_projection,
			       Coordinate destination_ul,
			       double destination_pixel_size);

	//! A constructor
	/*! 

	  This constructor takes six parameters. The first three
	  describe the source raster: a projection, upper-left
	  coordinate in projected coordinates, and a pixel size in meters.

	  The second set of three parameters are the same as the first
	  three but for the destination raster.

	*/
	RasterCoordTransformer(shared_ptr<Projection> source_projection,
			       Coordinate source_ul,
			       double source_pixel_size,
			       shared_ptr<Projection> destination_projection,
			       Coordinate destination_ul,
			       double destination_pixel_size);
	~RasterCoordTransformer();

	//! A normal member taking a single argument and returning an Area struct.
	/*

	  This function takes a coordinate in the source raster space
	  and maps is to an area in the destination raster space. The
	  returned area consists of two points: an upper-left
	  coordinate, and a lower-right. The coordinates are
	  _inclusive_. So if the area consists of a single point, the
	  upper-left == lower-right.

	   \param source a Coordinate struct that specifies the point in the source raster space to map to the destination raster space.
	*/
	Area Transform(Coordinate source);
private:
	shared_ptr<Projection> src_proj, dest_proj;
	Coordinate source_ul_;
	double source_pixel_size_;
	Coordinate destination_ul_;
	double destination_pixel_size_;
};

/**
 * This function creates a new raster file at the path
 * output_filename, with projection specified by output_srs. The
 * minbox of in is calculated with the new projection.
 *
 *
 * @param in The ProjectedRaster that is used to determine the size of the new raster
 * @param output_filename The path the new raster will be created at.
 * @param output_srs The Proj.4 specification of  projection and projection parameters.
 *
 */
PRB_ERROR CreateOutputRaster(shared_ptr<ProjectedRaster> in,
			string output_filename,
			double output_pixel_size,
			string output_srs);

/**
 * This function creates a new raster with a maximum size of
 * output_size x output_size, meant to be used to make reprojection previews.
 *
 * @param input The ProjectedRaster that is used to determine the geographic area  of the output file
 * @param output_filename File path that the new raster will be created at
 * @param output_pixel_size Size in meters of the pixels in the output raster
 * @param output_srs Projection specification string
 * @param Maximum pixel count of one dimension of new raster, eg 100 means raster will be maximum of 100x100
 *
 */
PRB_ERROR CreateSampleOutput(shared_ptr<ProjectedRaster> input,
			string output_filename,
			string output_srs, 
			int output_size);

/**
 * This function partitions the _area_ of the specified ProjectedRaster into approximately partition_count pieces.
 * 
 *
 * @param destination The ProjectedRaster to partition
 * @param partition_count The number of partitions to create
 * @return A std:vector of Areas. The Areas represent areas in raster coordinates. They will cover the entire
 *         raster and will not overlap. If the ul.x value of an Area is -1, the area is invalid and should be
 *         ignored. This will happen if you, for example, ask for 100 paritions from a raster with 99 pixels.  
 */
std::vector<Area> PartitionByCount(shared_ptr<ProjectedRaster> destination,
				   int partition_count);

Area ProjectedMinbox(shared_ptr<ProjectedRaster> input,
		     shared_ptr<Projection> output_projection);


Area RasterMinbox(shared_ptr<ProjectedRaster> source,
				shared_ptr<ProjectedRaster> destination,
				Area destination_raster_area);

bool ReprojectChunk(RasterChunk::RasterChunk *source, RasterChunk::RasterChunk *destination, string fillvalue, string resampler_name);

template <class pixelType>
bool ReprojectChunkType(RasterChunk::RasterChunk *source, RasterChunk::RasterChunk *destination, pixelType fillvalue, pixelType (*resampler)(Coordinate, 
																	     Coordinate,
																	     int,
																	     pixelType*))
{

	shared_ptr<Projection> outproj, inproj;
        Coordinate temp1, temp2;
	std::vector<char> inraster, outraster;

        outproj = destination->projection_;
        inproj = source->projection_;
	
	Area pixelArea;

	RasterCoordTransformer rt(destination->projection_, 
				  destination->ul_projected_corner_,
				  destination->pixel_size_,
				  source->projection_,
				  source->ul_projected_corner_,
				  source->pixel_size_);        


	for (int chunk_y = 0; chunk_y < destination->row_count_; ++chunk_y)  {
		for (int chunk_x = 0; chunk_x < destination->column_count_; ++chunk_x) {
			temp1.x = chunk_x; 
			temp1.y = chunk_y;
			
			pixelArea = rt.Transform(temp1);

			if (pixelArea.ul.x == -1.0) {
			  // The pixel is outside of the projected area
			  reinterpret_cast<pixelType*>(destination->pixels_)[chunk_x + chunk_y * destination->column_count_] = fillvalue;
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

			if (lr_x > (source->column_count_ - 1)) {
				lr_x = source->column_count_ - 1;
			}

			if (ul_y > (source->row_count_ - 1)) {
				ul_y = source->row_count_ - 1;
			}
			
			// Perform resampling...
			if (resampler != NULL && ((ul_x <= lr_x) || (lr_y <= ul_x))) { // ul/lr do not enclose an area, use NN
				reinterpret_cast<pixelType*>(destination->pixels_)[chunk_x + chunk_y * destination->column_count_] = 
				  reinterpret_cast<pixelType*>(source->pixels_)[ul_x + ul_y * source->column_count_];
				continue;
			}

			reinterpret_cast<pixelType*>(destination->pixels_)[chunk_x + chunk_y * destination->column_count_] = 
				resampler(Coordinate(ul_x, ul_y, UNDEF), 
							   Coordinate(lr_x, lr_y, UNDEF), 
							   source->column_count_,
							   reinterpret_cast<pixelType*>(source->pixels_));
		}
	}

	return true;
}

#endif // REPROJECTOR_HH
