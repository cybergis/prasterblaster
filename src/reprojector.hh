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

#ifndef REPROJECTOR_HH
#define REPROJECTOR_HH

#include <vector>

#include "gctp_cpp/projection.h"

#include "projectedraster.hh"
#include "resampler.hh"

using resampler::resampler_func;


class RasterCoordTransformer
{
public:
	RasterCoordTransformer(ProjectedRaster *source, ProjectedRaster *dest);
	~RasterCoordTransformer();
	Area Transform(Coordinate source);
private:
	ProjectedRaster *src, *dest;
	Projection *src_proj, *dest_proj;
};

class Reprojector
{
public:
	/* Constructor that creates a Reprojector object give an input and output raster.
	 *
	 * @param input Source raster
	 * @param output Target raster
	 */
	Reprojector(ProjectedRaster *_input, 
		    ProjectedRaster *_output);
	~Reprojector();

	/* 
	 * Initiates serial reprojection.
	 * 
	 * Only one processor is used.
	 */
	void reproject();

	/*
	 * Initiates parallel reprojection.
	 *
	 * The number and identity of nodes involved is determined by the mpi configuration.
	 */
	void parallelReproject();


	long startIndex(long process_number, vector<long> process_sizes);
	vector<long> getChunkSizes(long row_count, long chunk_count);
	vector<long> getChunkAssignments(long chunk_count, long process_count);
	void reprojectChunk(int firstRow, int numRows);
	int numprocs, rank;
	double maxx, minx, maxy, miny;
	ProjectedRaster *input;
	ProjectedRaster *output;

	resampler_func resampler;
};
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
