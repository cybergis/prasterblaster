// $Id: rasterinfo.h,v 1.5.2.3 2008/02/28 22:06:06 dmattli Exp $


#ifndef RASTERINFO_H
#define RASTERINFO_H

#include <QString>

class RasterInfo
{
public:
   //Constuctors and Destructor
   RasterInfo();
   RasterInfo( const QString &xmlFileName );
   RasterInfo( const RasterInfo &src );
   ~RasterInfo();

   //Files
   bool setXmlFileName( QString &xmlFileName );
   bool setImageFileName( QString &imageFilename );
   QString imgFileName() const;// {return fileName + ".img";}
   QString getXmlFileName() const;// {return fileName + ".xml";}

   //Author
   bool setAuthor( const QString &name, const QString &company, const QString &email );
   QString author() const {return aName;}
   QString company() const {return aCompany;}
   QString email() const {return aEmail;}

   //Area
   bool setArea( double ul_X, double ul_Y, int rows, int cols );
   bool setUL( double ul_X, double ul_Y );
   double ul_X() const {return ulx;}
   double ul_Y() const {return uly;}
   int rows() const {return row;}
   int cols() const {return col;}

   //Pixel Description
   bool setPixelDescription( const QString &dataType, double pixelSize, double fillValue, double noDataValue );
   bool setPixelDescription( bool isSigned, int bitsCount, const QString &type, double pixelSize, double fillValue, double noDataValue );
   bool setDataType( const QString &dataType );
   bool setDataType( bool isSigned, int bitsCount, const QString &type );
   bool setPixelSize( double pixelSize );
   bool setFillValue( double fillValue );
   bool setNoDataValue( double noDataValue );
   QString dataType() const {return datatype;}
   QString fullDataType() const;	//returns isSigned + bitCount + dataType 
   bool isSigned() const {return signd;}
   int bitCount() const {return bits;}
   QString type() const {return datatype;}
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
   bool load( QString xmlFileName );
   bool save( );
   bool save( QString xmlFileName );
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

   QString xmlFileName;
   QString imageFileName;

   QString  aName;
   QString  aCompany;
   QString  aEmail;

   //Used to keep hold of aName, aCompany, and aEmail between functions
   QString *tempAName, *tempACompany, *tempAEmail; 

   double   ulx;
   double   uly;
   int      row;
   int      col;

   bool     signd;

   int      bits;
   QString  datatype;
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
