// $Id: rasterinfo.cpp,v 1.9.2.3 2008/02/28 22:22:33 dmattli Exp $


#include <fstream>
#include <cstdio>
#include <sstream>

#include "rasterinfo.h"

#include "rasterxml.h"
#include "gctpnames.h"

using namespace std;

//Default constructor
// Nulls all private variables.
RasterInfo::RasterInfo()
{
     gctpParams = new double[15];

   defaults();
}

//Load constructor
// Loads the values found in the file whose name and path are
// given in imgFileName.
RasterInfo::RasterInfo( const string &xmlFileName )
{
     gctpParams = new double[15];

   if( !load( xmlFileName ) )
      defaults();
}

//Copy constructor
//
RasterInfo::RasterInfo( const RasterInfo &src )
{
     gctpParams = new double[15];

   copy( src );
}

//Destructor
// Deletes all the dynamic memory
RasterInfo::~RasterInfo()
{
         delete [] gctpParams;
}

bool RasterInfo::setXmlFileName( string &_xmlFileName )
{
  xmlFileName = _xmlFileName;

   return parseFileName();
}

bool RasterInfo::setImageFileName( string &_imageFileName )
{
  imageFileName = _imageFileName;

   return parseFileName();
}

string RasterInfo::imgFileName() const
{
   if( imageFileName == "" )
      return string();
   return imageFileName;
}

string RasterInfo::getXmlFileName() const
{
   if( xmlFileName == "" )
      return string();
   return xmlFileName;
}

bool RasterInfo::setAuthor( const string &name, const string &company, const string &email )
{
   aName = name.length()?name:"Unknown";
   aCompany = company.length()?company:"Unknown";
   aEmail = email.length()?email:"Unknown";

   tempAName = new string( aName );
   tempACompany = new string( aCompany );
   tempAEmail = new string( aEmail );

   return true;
}

bool RasterInfo::setArea( double ul_X, double ul_Y, int rows, int cols )
{
   ulx = ul_X;
   uly = ul_Y;
   row = rows;
   col = cols;

   return (row > 0 && col > 0);
}

bool RasterInfo::setUL( double ul_X, double ul_Y )
{
   ulx = ul_X;
   uly = ul_Y;

   return true;
}

bool RasterInfo::setPixelDescription( const string &dataType, double pixelSize, double fillValue, double noDataValue )
{
   return setDataType( dataType ) && setPixelSize( pixelSize ) && setFillValue( fillValue ) && setNoDataValue( noDataValue );
}

bool RasterInfo::setPixelDescription( bool isSigned, int bitsCount, const string &type, double pixelSize, double fillValue, double noDataValue )
{
   return setDataType( isSigned, bitsCount, type ) && setPixelSize( pixelSize ) && setFillValue( fillValue ) && setNoDataValue( noDataValue );
}

bool RasterInfo::setDataType( const string &dataType )
{
   //Handles "Unsigned/Signed 8/16/32 Bit Integer"
   //    and "Signed 32/64 Bit IEEE Float"

  if( dataType.find( "Signed" ) != string::npos )
      signd = true;
   else
      signd = false;

   if( dataType.find( "64" ) != string::npos )
      bits = 64;
   else if( dataType.find( "16" ) != string::npos )
      bits = 16;
   else if( dataType.find( "32" ) != string::npos )
      bits = 32;
   else //( dataType.contains( "8" ) > 0 )
      bits = 8;

   if( dataType.find("IEEE")  !=  string::npos
       || dataType.find("Float") !=  string::npos)
   {
      datatype = "IEEE Float";
      signd = true;
   }
   else
      datatype = "Integer";

   return !dataType.empty() && bits != 0;
}

bool RasterInfo::setDataType( bool isSigned, int bitsCount, const string &type )
{
   signd = isSigned;

   if( bitsCount % 8 == 0 )
      bits = bitsCount;
   else
      bits = 8;

   if( type == "Float" || type == "IEEE Float" || type == "IEEE" )
      datatype = "IEEE Float";
   else
      datatype = "Integer";

   return bitsCount > 0 && !type.empty();
}

bool RasterInfo::setPixelSize( double pixelSize )
{
   pixsize = pixelSize;

   return pixsize > 0;
}

bool RasterInfo::setFillValue( double fillValue )
{
   if( fillValue < 0 && !signd )
   {
      fillval = 0;
      return false;
   }

   fillval = fillValue;
   return true;
}

bool RasterInfo::setNoDataValue( double noDataValue )
{
   if( noDataValue < 0 && !signd )
   {
      noval = 0;
      return false;
   }

   noval = noDataValue;
   return true;
}


string RasterInfo::fullDataType() const
{
        string ret;
        std::stringstream out;
	
	ret = signd ? "Signed " : "Unsigned ";
        out << bits;
ret += out.str() + "Bit ";
	ret += datatype;

	return ret;	
}


bool RasterInfo::setHasFillValue( const bool& hasFill )
{
   hasFillVal = hasFill;

   return hasFillVal == hasFill;
}

bool RasterInfo::setHasNoDataValue( const bool& hasNoData )
{
   hasNoDataVal = hasNoData;

   return hasNoDataVal == hasNoData;
}

bool RasterInfo::setProjection( int projNumber, int zoneNumber, int datumNumber, int unitNumber )
{
   datumcode = datumNumber; unitcode = unitNumber;
   return setProjectionNumber( projNumber ) && setZoneNumber( zoneNumber );
}

bool RasterInfo::setProjectionNumber( int projNumber )
{
   if( projNumber >= 0 && projNumber < 31 )
      projcode = projNumber;

   return projcode == projNumber;
}

bool RasterInfo::setZoneNumber( int zoneNumber )
{
   zonecode = zoneNumber;

   if( ( projcode != 1 && projcode !=4 ) || //If not UTM or Lambert Conf. con.
      (zonecode > 60 && zonecode != 62) ||
      ( zonecode < -60 ) || ( zonecode == 0 ) )
      zonecode = 62;

   return zonecode == zoneNumber;
}

bool RasterInfo::setUnitNumber( int unitNumber )
{
   unitcode = unitNumber;



   return unitcode == unitNumber;
}

bool RasterInfo::setGctpParam(int i, double value)
{
   if( i >= 0 && i < 15 )
   {
      gctpParams[i] = value;
      return true;
   }

   return false;
}

bool RasterInfo::load()
{
   return load( xmlFileName );
}

bool RasterInfo::save()
{
   return save( xmlFileName );
}

bool RasterInfo::load( string _xmlFileName )
{

  ifstream ifile(_xmlFileName.c_str(), ifstream::in);

   defaults();
   xmlFileName = _xmlFileName;
   parseFileName();

   if(ifile)
   {
      return loadXml();
   }
   else
   {
      return false;
   }

   return true;
}

bool RasterInfo::save( string _xmlFileName )
{
   bool returnValue = false;

   if( !_xmlFileName.empty() )
   {
     xmlFileName = _xmlFileName;
     parseFileName();
   }

   try
   {
      RasterXML r;

	  r.setAuthorName( tempAName->c_str() );
	  r.setAuthorCompany( tempACompany->c_str() );
	  r.setAuthorEmail( tempAEmail->c_str() );

      r.setUlCorner( ulx, uly );
      r.setRows( row );
      r.setCols( col );

      r.setSigned( signd );
      r.setBits( bits );
      r.setDataType( datatype.c_str() );
      r.setPixelSize( pixsize );
      r.setFillValue( fillval );
      r.setNoDataValue( noval );

      r.setHasFillValue( hasFillVal );
      r.setHasNoDataValue( hasNoDataVal );

      r.setProjNumber( projcode );
      r.setProjName( projNames[projcode] );
      r.setZone( zonecode );
      r.setDatumNumber( datumcode );
      r.setDatumName( spheroidNames[datumcode] );
      r.setUnitsNumber( unitcode );
      r.setUnitsName( unitNames[unitcode] );

      char variation = 'a';
      if( projcode == 8 && gctpParams[8] == 1 )
         variation = 'b';
      if( ( projcode == 20 || projcode == 22 ) && gctpParams[12] == 1 )
         variation = 'b';
      vector<string> paramNames( gctpNames(projcode, variation) );
      for( int i = 0; i < 15; ++i )
         r.setParam( i, gctpParams[i], nameMeanings(paramNames[i]).c_str() );

      returnValue =  r.save( string( xmlFileName ).c_str() );
   }
   catch( ...  )
   {
     //      QMessageBox::critical( NULL, "Error", exception.getMessage() );
      returnValue = false;
   }

   return returnValue;
}


bool RasterInfo::parseFileName()
{
//	fileName.remove( ' ' ); // This is broken. Many real paths can have spaces!
/*
    if( fileName.right(4) == ".img" || fileName.right(4) == ".xml")
      fileName.truncate( fileName.length() - 4 );

   else if( fileName.right(9) == ".img.info" )
      fileName.truncate( fileName.length() - 9 );

   return !fileName.isNull();
*/
   return true;
}

void RasterInfo::loadInfo()
{
  /*
   QFile *file = new QFile( imageFileName + ".img.info" );
   file->open( QIODevice::ReadOnly );
   stringList inFile( string( file->readAll() ).split('\n') );
   file->close();
   delete file;

   ////////Rows and Columns
   int breakPoint = inFile[0].indexOf( ' ', 0, Qt::CaseInsensitive );
   row = inFile[0].left( breakPoint ).toInt();
   col = inFile[0].right( inFile[0].length() - breakPoint - 1 ).toInt();

   projcode = inFile[1].toInt(); //Projection Number/Name
   zonecode = inFile[2].toInt(); //Zone Code
   unitcode = 2; //Unit Type  NOTE: mapimg currently only supports meters
   datumcode = 19; //Spheroid Name  NOTE: mapimg currently only supports Sphere
   pixsize = inFile[5].toDouble(); //Pixel Size

   ////////UL Latitude and Longitude
   breakPoint = inFile[6].indexOf( ' ', 0, Qt::CaseInsensitive );
   ulx = inFile[6].left( breakPoint ).toDouble();
   uly = inFile[6].right( inFile[6].length() - breakPoint - 1 ).toDouble();

   ////////15 GCTP Params
   stringList gctpValues = string( inFile[7] ).split( " " );
   for( int i = 0; i < 15; ++i )
      gctpParams[i] = gctpValues[i].toDouble();

   setDataType( inFile[8] ); //Data Type
  */
}

bool RasterInfo::loadXml()
{
   bool returnValue = false;

   try
   {
      RasterXML xml( string( xmlFileName ).c_str() );

      aName = xml.getAuthorName();
      aCompany = xml.getAuthorCompany();
      aEmail = xml.getAuthorEmail();

      ulx = xml.getUlx();
      uly = xml.getUly();
      row = xml.getRows();
      col = xml.getCols();

      signd = xml.isSigned();
      bits = xml.getBits();
      datatype = xml.getDataType();
      pixsize = xml.getPixelSize();

      hasFillVal = xml.hasFillValue();
      fillval = xml.getFillValue();
      hasNoDataVal = xml.hasNoDataValue();
      noval = xml.getNoDataValue();

      projcode = xml.getProjNumber();
      zonecode = xml.getZone();
      datumcode = xml.getDatumNumber();
      unitcode = xml.getUnitsNumber();

      for( int i = 0; i < 15; ++i )
         gctpParams[i] = xml.getGCTPParam( i );

      returnValue = true;
   }
   catch( ... )
   {
     //      QMessageBox::critical( NULL, "XML File Error", exception.getMessage() );
      returnValue = false;
   }

   return returnValue;
}

bool RasterInfo::remove()
{
  ifstream file( xmlFileName.c_str() );
   if( file ) {
     file.close();
     return std::remove(xmlFileName.c_str());
   }

   return true;
}

bool RasterInfo::saveToTiff()
{

    return true;
}

bool RasterInfo::ready()
{
   if( !(row > 0 && col > 0) )
      return false;

   if( bits == 0 )
      return false;

   if( datatype.empty() )
      return false;

   if( !(pixsize > 0) )
      return false;

   if( projcode < 0 || projcode > 30 )
      return false;
   if( zonecode <= 0 || ( zonecode > 60 && zonecode != 62 ) )
      return false;

   if( datumcode != 19 || unitcode != 2 )
      return false;

   return true;
}

bool RasterInfo::notReady()
{
   return !ready();
}

void RasterInfo::defaults()
{
   xmlFileName = "";

   aName = "";  // "" author
   aCompany = "";
   aEmail = "";

   row = 0;                // 0x0
   col = 0;

   ulx = -1.0;             // (-1,-1)
   uly = -1.0;

   signd = false;          // Unsigned 8 Bit Integer
   bits = 8;
   datatype = "Integer";

   pixsize = 55597.454840; //30 Minutes

   hasFillVal = false;     // Undefined
   fillval = 0;

   hasNoDataVal = false;   // Undefined
   noval = 0;

   projcode = 17;          // Default to Equirectangular
   zonecode = 62;
   datumcode = 19;         // Sphere of Radius...
   unitcode = 2;           // Meters

   for( int i = 0; i < 15; ++i )
      gctpParams[i] = 0.0; // All 0's
}

void RasterInfo::copy( const RasterInfo &src )
{
   xmlFileName = src.xmlFileName;
   imageFileName = src.imageFileName;

   aName = src.aName;
   aCompany = src.aCompany;
   aEmail = src.aEmail;

   ulx = src.ulx;
   uly = src.uly;
   row = src.row;
   col = src.col;

   signd = src.signd;
   bits = src.bits;
   datatype = src.datatype;
   pixsize = src.pixsize;
   fillval = src.fillval;
   noval = src.noval;

   hasFillVal = src.hasFillVal;
   hasNoDataVal = src.hasNoDataVal;

   projcode = src.projcode;
   zonecode = src.zonecode;
   datumcode = src.datumcode;
   unitcode = src.unitcode;

   for( int i = 0; i < 15; ++i )
      gctpParams[i] = src.gctpParams[i];
}

bool RasterInfo::fakeIt()
{
   datumcode = 19; unitcode = 2;

   if( datatype.empty() )
   {
      signd = true;
      bits = 8;
      datatype = "Integer";
   }

   if( zonecode <= 0 || zonecode > 60 )
      zonecode = 62;

   return ready();
}
