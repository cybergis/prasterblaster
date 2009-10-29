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

	ProjectedRaster(string filename);
	ProjectedRaster(long num_rows, long num_cols, long pixel_bits,
			Projection *proj, double ulx, double uly); 
	~ProjectedRaster();


	void* getData();
	Projection* getProjection();
	ProjectedRaster* getSubRaster(int fromRow, int upToRow);
        bool isReady();
	bool write(string filename);
	
	// Area
	void setUL(double ul_x, double ul_y);
	void setRowCount(int rows);
	void setColCount(int cols);
	int getRowCount();
	int getColCount();
	
	// Pixel description
	void setDataType(const std::string &datatype);
	void setPixelSize(double pixsize);
	void setBitCount(int bits);
	void setSigned();
	void setUnsigned();
	void setInteger();
	void setFloat();
	string getDataType();
	double getPixelSize();
	int getBitCount();
	bool isSigned();
	bool isInteger();
	bool isFloat();
	
	
	// Projection
	void setProjection(ProjCode p);
	void setZoneNumber(int zone);
	void setUnit(ProjUnit unit);
	void setDatum(ProjDatum pd);
	void setGctpParams(double params[]);
	int getZoneNumber();
	ProjDatum getDatum();
	double* getGctpParams();

	// Members
	double ul_x, ul_y;

	// Pixel Description
	int rows, cols, bitCount;
	int start_row;
	bool issigned;
  	bool integer;
	double pixsize; // In degrees
	std::string type;
	std::string filename;

	// Projection
	int zoneNumber;
	ProjCode projectionCode;
	ProjDatum datum;
	ProjUnit unit;
	double gctpParams[16];
	
	void* data;
private:
	bool readImgRaster(string filename);
	bool readRaster(string filename);
	Projection *projection;
	Transformer t;
        ProjectedRaster& operator=(ProjectedRaster& a);
	bool ready;
};

#endif //RASTER_HH
