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

#include <boost/shared_ptr.hpp>

#include "gctp_cpp/projection.h"

#include "projectedraster.hh"
#include "resampler.hh"

using boost::shared_ptr;
using resampler::resampler_func;


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

class Chunker
{
public:
	Chunker(shared_ptr<ProjectedRaster> source, 
		shared_ptr<ProjectedRaster> destination);
	
	void clampGeoCoordinate(Coordinate *c);
	vector<ChunkExtent> getChunksByCount(int process_count, int process_index = -1);
	vector<ChunkExtent> getChunksBySize(int max_size, int process_index = -1);

	
private:
	shared_ptr<ProjectedRaster> source_raster;
	shared_ptr<ProjectedRaster> dest_raster;
		

};

class Reprojector
{
public:
	/* Constructor that creates a Reprojector object, give an input and output raster.
	 *
	 * @param input Source raster
	 * @param output Target raster
	 * @param number of process to launch
	 */ 
	Reprojector(shared_ptr<ProjectedRaster> _input, 
		    shared_ptr<ProjectedRaster> _output,
		    int  processCount,
		    int rank);
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

	void reprojectChunk(ChunkExtent chunk);
	int numprocs, rank, numchunks;
	double maxx, minx, maxy, miny;
	shared_ptr<ProjectedRaster> input;
	shared_ptr<ProjectedRaster> output;

	resampler_func resampler;
};

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
