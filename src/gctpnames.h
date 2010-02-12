// $Id: gctpnames.h,v 1.6.2.2 2008/01/15 19:58:42 dmattli Exp $


#ifndef GCTPNAMES_H
#define GCTPNAMES_H

#include <string>
#include <vector>

#include <QStringList>

using namespace std;

namespace
{
   //    gctpNames() returns a QStringList containing the names 
   // of each of the gctp parameters for the projection defined
   // by projNum and variation.
   vector<string> gctpNames( uint projNum, char variation )
   {
     //      QStringList retList;
      vector<string> retList;

      switch( projNum )
      {
      case 0:  //Geographic
	retList.push_back("NoWay");
         break;

      case 1:  //UTM
//         retList << "Lon/Z" << "Lat/Z" << "NoWay"; // GCTP doesn't need lat/lon anymore.
         break;

      case 2:  //State Plane
	retList.push_back("Sphere"); 
	retList.push_back("NoWay");
	break;

      case 3:  //Albers Equal Area
	retList.push_back("SMajor"); 
	retList.push_back("SMinor"); 
	retList.push_back("STDPR1"); 
	retList.push_back("STDPR2"); 
	retList.push_back("CentMer"); 
	retList.push_back("OriginLat"); 
	retList.push_back("FE"); 
	retList.push_back("FN"); 
         break;

      case 4:  //Lambert Conformal Conic
	retList.push_back("SMajor"); 
	retList.push_back("SMinor"); 
	retList.push_back("STDPR1"); 
	retList.push_back("STDPR2"); 
	retList.push_back("CentMer"); 
	retList.push_back("OriginLat"); 
	retList.push_back("FE"); 
	retList.push_back("FN"); 
         break;

      case 5:  //Mercator
	retList.push_back("SMajor"); 
	retList.push_back("SMinor"); 
	retList.push_back(""); 
	retList.push_back(""); 
	retList.push_back("CentMer"); 
	retList.push_back("TrueScale"); 
	retList.push_back("FE"); 
	retList.push_back("FN"); 
         break;

      case 6:  //PolarStereographic
	retList.push_back("SMajor"); 
	retList.push_back("SMinor"); 
	retList.push_back(""); 
	retList.push_back(""); 
	retList.push_back("LongPol"); 
	retList.push_back("TrueScale"); 
	retList.push_back("NoWay"); 
         break;

      case 7:  //Polyconic
	 	retList.push_back("SMajor"); 
	 	retList.push_back("SMinor"); 
	 	retList.push_back(""); 
	 	retList.push_back(""); 
	 	retList.push_back("CentMer"); 
	 	retList.push_back("OriginLat"); 
	 	retList.push_back("FE"); 
	 	retList.push_back("FN"); 
         break;

      case 8:  //Equid. Conic
	if( variation == 'a' ) {
	 retList.push_back("SMajor"); 
	 retList.push_back("SMinor"); 
	 retList.push_back("STDPAR"); 
	 retList.push_back(""); 
	 retList.push_back("CentMer"); 
	 retList.push_back("OriginLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("zero"); 
	} else {
	   retList.push_back("SMajor"); 
	 retList.push_back("SMinor"); 
	 retList.push_back("STDPAR"); 
	 retList.push_back(""); 
	 retList.push_back("CentMer"); 
	 retList.push_back("OriginLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("zero"); 
	}
         break;

      case 9: //Transverse Mercator
	 retList.push_back("SMajor"); 
	 retList.push_back("SMinor"); 
	 retList.push_back("FactorM"); 
	 retList.push_back(""); 
	 retList.push_back("CentMer"); 
	 retList.push_back("OriginLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
         break;

      case 10: //Stereographic
	 retList.push_back("Sphere"); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back("CentLon"); 
	 retList.push_back("CenterLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("NoWay"); 
         break;

      case 11: //Lambert Azimuthal
      case 12: //Azimuthal
      case 14: //Orthographic
	 retList.push_back("Sphere"); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back("CentLon"); 
	 retList.push_back("CenterLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
         break;

      case 13: //Gnomonic
	 retList.push_back("Sphere"); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back("CentLon"); 
	 retList.push_back("CenterLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("NoWay"); 
         break;

      case 15: //Gen. Vert. Near Per
	 retList.push_back("Sphere"); 
	 retList.push_back(""); 
	 retList.push_back("Height"); 
	 retList.push_back(""); 
	 retList.push_back("CentLon"); 
	 retList.push_back("CenterLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
         break;

      case 16: //Sinusoidal
      case 18: //Miller Cylindrical
      case 21: //Robinson
      case 25: //Mollweide
      case 27: //Hammer
      case 28: //Wagner IV
      case 29: //Wagner VII
	 retList.push_back("Sphere"); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back("CentMer"); 
	 retList.push_back(""); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
         break;

      case 17: //Equirectangular
	 retList.push_back("Sphere"); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back("CentMer"); 
	 retList.push_back("TrueScale"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
         break;

      case 19: //Van der Grinten
	 retList.push_back("Sphere"); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back("CentMer"); 
	 retList.push_back("OriginLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
         break;

      case 20: //Hotin Oblique Merc
	if( variation == 'a' ) {
	 retList.push_back("SMajor"); 
	 retList.push_back("SMinor"); 
	 retList.push_back("FactorH"); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back("OriginLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("Long1"); 
	 retList.push_back("Lat1"); 
	 retList.push_back("Long2"); 
	 retList.push_back("Lat2"); 
	 retList.push_back("zero"); 
	} else {
	 retList.push_back("SMajor"); 
	 retList.push_back("SMinor"); 
	 retList.push_back("FactorH"); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back("OriginLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("Long1"); 
	 retList.push_back("Lat1"); 
	 retList.push_back("Long2"); 
	 retList.push_back("Lat2"); 
	 retList.push_back("zero"); 
	}
         break;

      case 22: //Space Oblique Merc
	if( variation == 'a' ) {
	 retList.push_back("SMajor"); 
	 retList.push_back("SMinor"); 
	 retList.push_back(""); 
	 retList.push_back("IncAng"); 
	 retList.push_back("AscLong"); 
	 retList.push_back(""); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("PSRev"); 
	 retList.push_back("LRat"); 
	 retList.push_back("PFlag"); 
	 retList.push_back(""); 
	 retList.push_back("zero"); 
	} else {
 	 retList.push_back("SMajor"); 
	 retList.push_back("SMinor"); 
	 retList.push_back(""); 
	 retList.push_back("IncAng"); 
	 retList.push_back("AscLong"); 
	 retList.push_back(""); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("PSRev"); 
	 retList.push_back("LRat"); 
	 retList.push_back("PFlag"); 
	 retList.push_back(""); 
	 retList.push_back("one"); 
	}
         break;

      case 23: //Alaska Conformal
	 retList.push_back("SMajor"); 
	 retList.push_back("SMinor"); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back(""); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("NoWay"); 
         break;

      case 24: //Interrupted Goode
      case 26: //Interrupted Mollweide

         break;

      case 30: //Oblated Equal Area
	 retList.push_back("Sphere"); 
	 retList.push_back(""); 
	 retList.push_back("Shapem"); 
	 retList.push_back("Shapen"); 
	 retList.push_back("CentLon"); 
	 retList.push_back("CenterLat"); 
	 retList.push_back("FE"); 
	 retList.push_back("FN"); 
	 retList.push_back("Angle"); 
         break;
      }

      //Fill the list up to 15.
      while( retList.size() < 15 )
	retList.push_back("");

      return retList;
   }

   // Extended names with more meaning then their gctpNames code
   string nameMeanings( string &gctpName )
   {
      if( gctpName == "" )
         return "Unused";
      else if( gctpName == "Angle" )
         return "Oval Rotation Angle:";
      else if( gctpName == "AscLong" )	
         return "Longitude of Ascending Orbit:";
      else if( gctpName == "AziAng" )	
         return "Azimuth Angle East of North:";
      else if( gctpName == "AzmthPt" )		
         return "Longitude where Azimuth Occurs:";
      else if( gctpName == "CenterLat" )	
         return "Latitude of Projection Center:";
      else if( gctpName == "CentLon" )		
         return "Longitude of Projection Center:";
      else if( gctpName == "CentMer" )	
         return "Longitude of Central Meridian:";
      else if( gctpName == "FactorH" )	
         return "Scale Factor at Center of Projection:";
      else if( gctpName == "FactorM" )		
         return "Scale Factor at Central Meridian:";
      else if( gctpName == "FE" )		
         return "False Easting:";
      else if( gctpName == "FN" )	
         return "False Northing:";
      else if( gctpName == "Height" )	
         return "Height of Perspective:";
      else if( gctpName == "IncAng" )	
         return "Inclination of Orbit:";
      else if( gctpName == "Lat/Z" )	
         return "Latitude within Zone:";
      else if( gctpName == "Lat1" )		
         return "Latitude of 1st Point on Center Line:";
      else if( gctpName == "Lat2" )
         return "Latitude of 2nd Point on Center Line:";
      else if( gctpName == "Lon/Z" )	
         return "Longitude within Zone:";
      else if( gctpName == "Long1" )	
         return "Longitude of 1st Point on Center Line:";
      else if( gctpName == "Long2" )	
         return "Longitude of 2nd Point on Center Line:";
      else if( gctpName == "LongPol" )	
         return "Longitude down below Pole of Map:";
      else if( gctpName == "LRat" )	
         return "Landsat Ratio:";
      else if( gctpName == "OriginLat" )	
         return "Latitude of Projection Origin:";
      else if( gctpName == "Path" )		
         return "Landsat Path Number:";
      else if( gctpName == "PFlag" )	
         return "End of Path Flag:";
      else if( gctpName == "PSRev" )	
         return "Period of Satellite Revolution:";
      else if( gctpName == "Satnum" )	
         return "Landsat Satellite Number:";
      else if( gctpName == "Shapem" )		
         return "Oval Shape Parameter m:";
      else if( gctpName == "Shapen" )	
         return "Oval Shape Parameter n:";
      else if( gctpName == "SMajor" )		
         return "Semi-Major Axis:";
      else if( gctpName == "SMinor" )	
         return "Semi-Minor Axis:";
      else if( gctpName == "Sphere" )	
         return "Radius of the Reference Sphere:";
      else if( gctpName == "STDPAR" )	
         return "Latitude of Standard Parallel:";
      else if( gctpName == "STDPR1" )	
         return "Latitude of 1st Standard Parallel:";
      else if( gctpName == "STDPR2" )	
         return "Latitude of 2nd Standard Parallel:";
      else if( gctpName == "TrueScale" )		
         return "Latitude of True Scale:";
      else if( gctpName == "one" )	
         return "One";
      else if( gctpName == "zero" )
         return "Zero";

      return "Unused";
   }


   // Convert the index of the projection combo box to actual projection code
   uint combo2proj( uint i )
   {
      if( i < 9 ) return i - 1;
	   else if( i < 20 ) return i;
	   else if( i == 20 ) return i + 1;
	   else if( i < 29 ) return i + 2;
	   else if( ( i == 29 ) || ( i == 30 ) ) return 22;
	   else if( ( i == 31 ) || ( i == 32 ) ) return 20;
	   return 8;
   }

   // A list of projections sorted by projection code
   const char *projNames[] =
     {
       "Geographic,Universal Transverse Mercator,State Plane Coordinates",
       "Albers Conical Equal Area",
       "Lambert Conformal Conic", 
       "Mercator",
       "Polar Stereographic",
       "Polyconic,Equidistant Conic", 
       "Transverse Mercator",
       "Stereographic",
       "Lambert Azimuthal Equal Area",
       "Azimuthal Equidistant", 
       "Gnomonic",
       "Orthographic",
       "General Vertical Near-Side Perspective",
       "Sinusoidal,"
       "Equirectangular",
       "Miller Cylindrical", "Van der Grinten",
       "Hotine Oblique Mercator", 
       "Robinson",
       "Space Oblique Mercator",
       "Modified Stereographic Conformal--Alaska",
       "Interrupted Goode Homolsine",
       "Mollweide",
       "Interrupted Mollweide", 
       "Hammer", 
       "Wagner IV",
       "Wagner VII",
       "Oblated Equal Area" 
     };


   // A list of unit types sorted by unit-code
   const char *unitNames[] = 
     {
       "Radians",
       "U.S. Feet",
       "Meters", 
       "Seconds of Arc", 
       "Degrees of Arc",
       "International Feet",
       "State Zone Table",
       "Unknown"
     };

   // A list of Spheroid types sorted by spheroid-code
   const char *spheroidNames[] = 
     {
       "Clarke 1866",
       "Clarke 1880",
       "Bessel", 
       "International 1967", 
       "International 1909",
       "WGS 72", 
       "Everest", 
       "WGS 66", 
       "GRS 1980", 
       "Airy", 
       "Modified Everest", 
       "Modified Airy",
       "WGS 84", 
       "Southeast Asia", 
       "Australlian National", 
       "Krassovsky", 
       "Hough", 
       "Mercury 1960",
       "Modified Mercury 1968", 
       "Sphere of Radius 6370997 meters", 
       "Unknown" 
     };

   // A list of common named pixel sizes sorted by decreasing values
   const char *pixelSizes[] =
     {
       "5 Degrees", 
       "1 Degree", 
       "30 Minutes", 
       "5 Minutes", 
       "30 Arc Seconds", 
       "Meters..."
     };

   // A list of the pixel values associated with the names
   const char *pixelValues[] = 
     {
       "555974.548395",
       "111194.909679", 
       "55597.454840", 
       "9266.242473", 
       "926.624247" 
     };

   // A list of all supported data types for mapimg
   const char *dataTypes[] = 
     {
       "Signed 64 Bit IEEE Float",
       "Signed 32 Bit IEEE Float",
       "Unsigned 32 Bit Integer",
       "Signed 32 Bit Integer",
       "Unsigned 16 Bit Integer", 
       "Signed 16 Bit Integer",
       "Unsigned 8 Bit Integer", 
       "Signed 8 Bit Integer" 
     };
}

#endif
