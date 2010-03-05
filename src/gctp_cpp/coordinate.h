#ifndef COORDINATE_H
#define COORDINATE_H

#include <math.h>
#include "constants.h"

/*! Coordinate struct
   This class provides a more readable way of storing and passing coordinate
   parameters for the Transformer class. It stores x and y as doubles and units
   corresponds to the GCTP enumeration as defined in constants.h
   */
struct Coordinate {
   /*! Default Constructor
      0's all attributes.
      */
   Coordinate();

   /*! Full Constructor
      Set's all attributes in the Coordinate according to the parameters.
      */
   Coordinate( double xx, double yy, ProjUnit uunits );

   /*! Copy Constructor
      Set's all attributes to equal those in Coordinate c.
      */
   Coordinate( const Coordinate &c );

   //! Set this Coordinate's attributes to equal those in Coordinate c.
   void copy( const Coordinate &c );

   double x;
   double y;
   ProjUnit units;
};

/*! Checks if two Coordinates are "close" to each other.
   Will return true if they have the same units and both x and y are
   within delta of each other.
   \param c1 First coordinate to compare
   \param c2 Second coordinate to compare
   \param delta Maximum difference allowed before the two coordinates are incomparable
   */
inline bool comparable( const Coordinate &c1, const Coordinate &c2, double delta )
{
   if( c1.units != c2.units ) return false;
   if( fabs(c1.x-c2.x) > delta ) return false;
   if( fabs(c1.y-c2.y) > delta ) return false;

   return true;
}

#endif
