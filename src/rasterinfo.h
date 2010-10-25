// $Id: rasterinfo.h,v 1.5.2.3 2008/02/28 22:06:06 dmattli Exp $


#ifndef RASTERINFO_H
#define RASTERINFO_H

#include <string>

using namespace std;

/*! ProjectedRaster class.
 *
 * This class parses and writes raster configuration files
 */

class RasterInfo
{
public:
   //Constuctors and Destructor
  //! Constructor
  /*!
   * This constructor takes no arguments. 
   */
   RasterInfo();

   //! Constructor
   /*!
    * This constructor takes a filename argument. The file specified
    * will be parsed.
    */ 
   RasterInfo( const string &xmlFileName );

   //! Constructor
   /*!
    * This constructor takes another RasterInfo object as an argument
    * and produces a copy.
    */
   RasterInfo( const RasterInfo &src );
   ~RasterInfo();

   //Files
   //! A normal member taking a string argument
   /*!
    * This function sets the filename of the associated xml configuration
    */
   bool setXmlFileName( string &xmlFileName );
   bool setImageFileName( string &imageFilename );
   string imgFileName() const;// {return fileName + ".img";}
   string getXmlFileName() const;// {return fileName + ".xml";}

   //Author
   bool setAuthor( const string &name, const string &company, const string &email );
   string author() const {return aName;}
   string company() const {return aCompany;}
   string email() const {return aEmail;}

   //Area
   bool setArea( double ul_X, double ul_Y, int rows, int cols );
   bool setUL( double ul_X, double ul_Y );
   double ul_X() const {return ulx;}
   double ul_Y() const {return uly;}
   int rows() const {return row;}
   int cols() const {return col;}

   //Pixel Description
   bool setPixelDescription( const string &dataType, double pixelSize, double fillValue, double noDataValue );
   bool setPixelDescription( bool isSigned, int bitsCount, const string &type, double pixelSize, double fillValue, double noDataValue );
   bool setDataType( const string &dataType );
   bool setDataType( bool isSigned, int bitsCount, const string &type );
   bool setPixelSize( double pixelSize );
   bool setFillValue( double fillValue );
   bool setNoDataValue( double noDataValue );
   string dataType() const {return datatype;}
   string fullDataType() const;	//returns isSigned + bitCount + dataType 
   bool isSigned() const {return signd;}
   int bitCount() const {return bits;}
   string type() const {return datatype;}
   double pixelSize() const {return pixsize;}
   double fillValue() const {return fillval;}
   double noDataValue() const {return noval;}

   bool hasFillValue() const {return hasFillVal;}
   bool hasNoDataValue() const {return hasNoDataVal;}
   bool setHasFillValue( const bool& hasFill );
   bool setHasNoDataValue( const bool& hasNoData );

   //Projection
   bool setProjection( int projNumber, int zoneNumber = 62, int datumNumber = 19, int unitNumber = 2 );
   bool setProjectionNumber( int projNumber );
   bool setZoneNumber( int zoneNumber );
   bool setUnitNumber( int unitNumber );
   int projectionNumber() const {return projcode;}
   int zoneNumber() const {return zonecode;}
   int datumNumber() const {return datumcode;}
   int unitNumber() const {return unitcode;}

   //GCTP Parameters
   bool setGctpParam(int i, double value);
   double gctpParam(int i) const {return ((i>=0)&&(i<15))?gctpParams[i]:0.0;}
    double *allGctpParams() const {return gctpParams;}

   //I/O
   bool load();
   bool load( string xmlFileName );
   bool save( );
   bool save( string xmlFileName );
   bool remove();
   bool saveToTiff();

   //Checks
   bool fakeIt();
   bool ready();
   bool notReady();

   void copy( const RasterInfo &src );


   void defaults();
   bool parseFileName();
   void loadInfo();
   bool loadXml();

   string xmlFileName;
   string imageFileName;

   string  aName;
   string  aCompany;
   string  aEmail;

   //Used to keep hold of aName, aCompany, and aEmail between functions
   string *tempAName, *tempACompany, *tempAEmail; 

   double   ulx;
   double   uly;
   int      row;
   int      col;

   bool     signd;

   int      bits;
   string  datatype;
   double   pixsize;
   double   fillval;
   double   noval;

   bool     hasFillVal;
   bool     hasNoDataVal;

   long     projcode;
   long     zonecode;
   long     datumcode;
   long     unitcode;

   double*  gctpParams;   
};

#endif //RASTERINFO_H
