/**
 * @file
 * @author  David Mattli
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * The ProjectedRaster class represents a raster with a location and a projection.
 */

#ifndef RASTER_HH
#define RASTER_HH

#include <string>

#include <gdal/gdal.h>
#include <gdal/gdal_priv.h>

#include "gctp_cpp/constants.h"
#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"

#include "rasterinfo.h"

using namespace std; // Don't do this :(

//! ProjectedRaster class.
/*!
  This class represents a raster with a projection and location.
*/

class ProjectedRaster
{
public:
	//! Constructor
/*! 
This constructor takes a single arguments, filename, representing
	the path to the raster to be opened.
*/

	ProjectedRaster(string _filename);

/*! 
  This constructor takes lots of arguments.
 */
	ProjectedRaster(string filename, 
			int num_rows, int num_cols,
			GDALDataType pixel_type, double pixel_size,
			int band_count,
			Projection *proj,
			double ulx, double uly);

	ProjectedRaster(string filename,
			ProjectedRaster *input,
			Projection *output_proj,
			GDALDataType pixel_type,
			double pixel_size);

/*!
  Destructor
 */
	~ProjectedRaster();


	Projection* getProjection();
        bool isReady();
	bool write(string filename);
	
	// Area
	int getRowCount();
	int getColCount();
	
	// Pixel description
	GDALDataType getPixelType();
	int bitsPerPixel();
	int bandCount();
	double getPixelSize();

	// Projection
	int getZoneNumber();
	ProjDatum getDatum();
	double* getGctpParams();

	// IO
	bool readRaster(int firstRow, int numRows, void* data);
	bool writeRaster(int firstRow, int numRows, void* data);

	// Members
	double ul_x, ul_y;
	int rows, cols;
	GDALDataType type;
	double pixel_size;
	int band_count;

	// File Description
	std::string filename;

	// Projection
	int zoneNumber;
	ProjCode projectionCode;
	ProjUnit unit;
	double gctpParams[16];
	
	void* data;
private: 
	bool loadImgRaster(string filename);
	bool readImgRaster(string filename);
	bool readRaster(string filename);
	Projection *projection;
	GDALDataset *dataset;
	Transformer t;
        ProjectedRaster& operator=(ProjectedRaster& a);
	bool ready;
};

#endif //RASTER_HH
