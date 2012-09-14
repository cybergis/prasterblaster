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

#ifndef SRC_PROJECTEDRASTER_H_
#define SRC_PROJECTEDRASTER_H_

#include <string>

#include <gdal.h>
#include <gdal_priv.h>

#include "gctp_cpp/coordinate.h"
#include "gctp_cpp/constants.h"
#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"

#include "rasterchunk.h"
#include "reprojection_tools.h"
#include "sharedptr.h"
#include "utils.h"

using std::string;

namespace librasterblaster  {

/*! ProjectedRaster class.
 *
 * This class represents a raster with a projection and location.
 */
class ProjectedRaster {
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
  shared_ptr<Projection> projection();

  //! A normal member taking no arguments.
  /*!
   * \return Returns true if the raster is in a good state for
   * reading/writing. If it returns false something is wrong and the
   * raster can't be trusted.
   */
  bool ready();

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
  int row_count();

  //! A normal member function taking no arguments.
  /*!
    \return Returns the geographical minbox of the raster.
  */
  Area geographical_minbox();
	
  //! A normal member function taking no arguments.
  /*!
    \return Returns the projected minbox of the raster.
  */
  Area projected_minbox();
	
  //! A normal member function taking no arguments.
  /*!
   * \return Returns the number of columns in the raster.
   */
  int column_count();
	
  // Pixel description
  //! A normal member function taking no arguments.
  /*!
   * \return Returns the datatype of the pixels
   */
  GDALDataType pixel_type();

  //! A normal member function taking no arguments.
  /*!
   * \return Returns the number of bits in each pixel
   */
  int bits_per_pixel();
  //! A normal member function taking no arguments.
  /*! 
   * \return Returns the number of bands in the raster
   */
  int band_count();
  //! A normal member function taking no arguments.
  /*!
   * \return Returns the size of the pixels in meters
   */
  double pixel_size();
 
  // Projection
  //! A normal member function taking no arguments.
  /*!
   * \return Returns the UTM zone number. If the projection is not UTM, this value is undefined.
   */
  int zone_number();
  //! A normal member function taking no arguments.
  /*!
   * \return Returns the ProjDatum enum value indicated the datum used for the raster's projection.
   */
  ProjDatum datum();

  //! A normal member function taking no arguments.
  /*!
   * \returns Returns a pointer to an array of the GCTP parameters
   */
  double* gctp_parameters();

  // IO
  //! A normal member function taking three arguments
  /*!
   * @param firstRow The index of the first row to be read
   * @param numRows The count of rows to be read
   * @param data Pointer to area to copy raster section
   * \returns A bool indicated a success or failure
   */
  bool read_raster(int firstRow, int numRows, void* data);

  //! A normal member function taking three arguments
  /*!
   * @param firstRow The index of the fist row to be written
   * @param numRows The count of rows to be written
   * @param data Pointer to data to be written
   * \returns A bool indicated a success or failure
   */
  bool write_raster(int firstRow, int numRows, void* data);

  //! A normal member function taking one argument
  /*!
   * @param area This represents the rectangle in raster coordinates
   * to be represented by the RasterChunk 
   * \returns A pointer to a
   *          RasterChunk, where the pixels are read from the ProjectedRaster
   */
  RasterChunk* create_raster_chunk(Area area);
  //! A normal member function taking one argument
  /*!
   * @param area This represents the rectangle in raster coordinates to be represented by the RasterChunk
   * \returns A pointer to a RasterChunk, where the pixels are not read, 
   *          but memory of a sufficient size is allocated
   */

  RasterChunk* create_allocated_raster_chunk(Area area);
  //! A normal member function taking one argument
  /*!
   * @param area This represents the rectangle in raster coordinates to be represented by the RasterChunk
   * \returns A pointer to a RasterChunk, where the pixels are not read, 
   *          and no memory is allocated.
   */

  RasterChunk* create_empty_raster_chunk(Area area);

  bool write_raster_chunk(RasterChunk *chunk);

  double ul_x();
  double ul_y();
  string srs();
 private:
  // Members
  double ul_x_, ul_y_;
  int row_count_, column_count_;
  GDALDataType pixel_type_;
  double pixel_size_;
  int band_count_;

  // File Description
  std::string filename_;

  // Projection
  int zone_number_;
  ProjCode projection_code_;
  ProjUnit unit_;
  double gctp_parameters_[16];
	
  void clamp_geo_coordinate(Coordinate *c);
  bool load_raster(string filename);
  static bool make_raster(string filename,
                         int cols,
                         int rows,
                         int band_count,
                         double ul_x,
                         double ul_y,
                         GDALDataType type,
                         shared_ptr<Projection> projection,
                         double pixel_size);

  shared_ptr<Projection> projection_;
  GDALDataset *dataset_;
  Area geographical_minbox_, projected_minbox_;
  Transformer t_;
  ProjectedRaster& operator=(ProjectedRaster& a);
  bool ready_;
};
}

#endif //  SRC_PROJECTEDRASTER_H_

