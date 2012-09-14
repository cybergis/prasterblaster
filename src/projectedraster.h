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
 * The ProjectedRaster class represents a raster with a location and a projection.
 *
 */

#ifndef RASTER_HH
#define RASTER_HH

#include <string>

#include <gdal.h>
#include <gdal_priv.h>

#include "gctp_cpp/coordinate.h"
#include "gctp_cpp/constants.h"
#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"

#include "rasterchunk.hh"
#include "sharedptr.hh"

using std::string;

struct Area {
	Area() {}
	Area(double ulx, double uly, double lrx, double lry) : ul(ulx, uly, UNDEF), lr(lrx, lry, UNDEF) {}
		
	Coordinate ul;
	Coordinate lr;
	ProjUnit units;
	
};

/*! ProjectedRaster class.
 *
 * This class represents a raster with a projection and location.
 */
class ProjectedRaster
{
public:
	//! Constructor
	/*! 
	 * This constructor takes a single argument, filename, representing
	 * the path to the raster to be opened.
	 */

	ProjectedRaster(string _filename);
	
	//! Constructor
	/*! 
	 * This constructor builds a ProjectedRaster from the arguments given.
	 * 
	 * @param filename The filename of the ProjectedRaster
	 * @param num_rows Number of rows
	 * @param num_cols Number of columns
	 * @param pixel_type The type of pixels
	 * @param pixel_size Size in meters of one of the pixel dimensions
	 * @param band_count Number of band in the raster
	 * @param proj Pointer to Projection object that describes the raster's projection
	 */
	static bool CreateRaster(string filename, 
				 int num_rows, int num_cols,
				 GDALDataType pixel_type, double pixel_size,
				 int band_count,
				 shared_ptr<Projection> proj,
				 double ulx, double uly);
	
	//! Constructor
	/*!
	 * This constructor creates a raster from a filename and an xml description file.
	 */
	static bool CreateRaster(shared_ptr<ProjectedRaster> input,
				 string filename,
				 string xmlDescriptionPath);
	
	static bool CreateRaster(string filename,
				 shared_ptr<ProjectedRaster> input,
				 shared_ptr<Projection> output_proj,
				 GDALDataType pixel_type,
				 double pixel_size);
	
 /*!
 * Destructor
 */
	~ProjectedRaster();


	//! A normal member taking no arguments.
/*!
 * \return A copy of the ProjectedRaster's projection object. It's the
 * callers responsibility to delete.
 */
  shared_ptr<Projection> getProjection();

//! A normal member taking no arguments.
/*!
 * \return Returns true if the raster is in a good state for
 * reading/writing. If it returns false something is wrong and the
 * raster can't be trusted.
 */
        bool isReady();

//! A normal member taking a single string argument
/*!
 * \param filename is a string indicating where to write the raster to.
 */
	bool write(string filename);
	
	/////// Area

	//! A normal member function taking no arguments.
/*! 
 * \return Returns the number of rows in the raster.
 *
 */
	int getRowCount();

  //! A normal member function taking no arguments.
  /*!
    \return Returns the geographical minbox of the raster.
  */
	Area getGeographicalMinbox();
	
	//! A normal member function taking no arguments.
  /*!
    \return Returns the projected minbox of the raster.
  */
	Area getProjectedMinbox();
	

	//! A normal member function taking two Arguments
	Coordinate getProjectedCoordinate(int rasterX, int RasterY);
	Coordinate getGeographicalCoordinate(int rasterX, int rasterY);

	int getPixelIndex(double longitude, double latitude);
	int getPixelIndex(Coordinate geographicalCoordinate);
  
  //! A normal member function taking no arguments.
/*!
 * \return Returns the number of columns in the raster.
 */
	int getColCount();
	
  // Pixel description
//! A normal member function taking no arguments.
/*!
 * \return Returns the datatype of the pixels
 */
	GDALDataType getPixelType();

//! A normal member function taking no arguments.
/*!
 * \return Returns the number of bits in each pixel
 */
	int bitsPerPixel();
//! A normal member function taking no arguments.
/*! 
 * \return Returns the number of bands in the raster
 */
	int bandCount();
//! A normal member function taking no arguments.
/*!
 * \return Returns the size of the pixels in meters
 */
	double getPixelSize();
 
	// Projection
//! A normal member function taking no arguments.
/*!
 * \return Returns the UTM zone number. If the projection is not UTM, this value is undefined.
 */
	int getZoneNumber();
//! A normal member function taking no arguments.
/*!
 * \return Returns the ProjDatum enum value indicated the datum used for the raster's projection.
 */
	ProjDatum getDatum();

//! A normal member function taking no arguments.
/*!
 * \returns Returns a pointer to an array of the GCTP parameters
 */
	double* getGctpParams();

  // IO
//! A normal member function taking three arguments
/*!
 * @param firstRow The index of the first row to be read
 * @param numRows The count of rows to be read
 * @param data Pointer to area to copy raster section
 * \returns A bool indicated a success or failure
 */
	bool readRaster(int firstRow, int numRows, void* data);

//! A normal member function taking three arguments
/*!
 * @param firstRow The index of the fist row to be written
 * @param numRows The count of rows to be written
 * @param data Pointer to data to be written
 * \returns A bool indicated a success or failure
 */
	bool writeRaster(int firstRow, int numRows, void* data);

  //! A normal member function taking one argument
  /*!
   * @param area This represents the rectangle in raster coordinates
   * to be represented by the RasterChunk 
   * \returns A pointer to a
   *          RasterChunk, where the pixels are read from the ProjectedRaster
   */
	RasterChunk::RasterChunk* createRasterChunk(Area area);
  //! A normal member function taking one argument
  /*!
   * @param area This represents the rectangle in raster coordinates to be represented by the RasterChunk
   * \returns A pointer to a RasterChunk, where the pixels are not read, 
   *          but memory of a sufficient size is allocated
   */

	RasterChunk::RasterChunk* createAllocatedRasterChunk(Area area);
  //! A normal member function taking one argument
  /*!
   * @param area This represents the rectangle in raster coordinates to be represented by the RasterChunk
   * \returns A pointer to a RasterChunk, where the pixels are not read, 
   *          and no memory is allocated.
   */

	RasterChunk::RasterChunk* createEmptyRasterChunk(Area area);
	bool writeRasterChunk(RasterChunk::RasterChunk *chunk);

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
	void clampGeoCoordinate(Coordinate *c);
	bool configureFromXml(string filename);
	bool loadImgRaster(string rasterFilename, string xmlFilename);
	bool loadRaster(string filename);
	static bool makeRaster(string filename,
			       int cols,
			       int rows,
			       int band_count,
			       double ul_x,
			       double ul_y,
			       GDALDataType type,
			       shared_ptr<Projection> projection,
			       double pixel_size);
	shared_ptr<Projection> projection;
	GDALDataset *dataset;
	Area geographicalMinbox, projectedMinbox;
	Transformer t;
        ProjectedRaster& operator=(ProjectedRaster& a);
	bool ready;
};


#endif //RASTER_HH
