

#ifndef CONST_H

#define CONST_H



#define PI 	3.141592653589793238

#define HALF_PI (PI*0.5)

#define TWO_PI 	(PI*2.0)

#define EPSLN	1.0e-10

#define R2D     57.2957795131



#define D2R     1.745329251994328e-2

#define S2R	4.848136811095359e-6



#define OK	0

#define ERROR  -1

#define IN_BREAK -2

#define MAX_VAL 4

#define MAXLONG 2147483647.

#define DBLLONG 4.61168601e18

//define NULL 0

#define RADVAL 19



#define LANDSAT_RATIO 0.5201613



/* The STPLN_TABLE unit value is specifically used for State Plane -- if units

   equals STPLN_TABLE and Datum is NAD83--actual units are retrieved from

   a table according to the zone.  If Datum is NAD27--actual units will be feet.

   An error will occur with this unit if the projection is not State Plane.  */



#define STPLN_TABLE 6



/* General code numbers */



#define IN_BREAK -2		/*  Return status if the interupted projection

				    point lies in the break area */

#define COEFCT 15		/*  projection coefficient count */

#define PROJCT 30		/*  projection count */

#define SPHDCT 31		/*  spheroid count */



#define MAXPROJ 31		/*  Maximum projection number */

#define MAXUNIT 5		/*  Maximum unit code number */

#define GEO_TERM 0		/*  Array index for print-to-term flag */

#define GEO_FILE 1		/*  Array index for print-to-file flag */

#define GEO_TRUE 1		/*  True value for geometric true/false flags */

#define GEO_FALSE -1		/*  False val for geometric true/false flags */

#define SQUARE(x)       x * x   /* x**2 */

#define CUBE(x)     x * x * x   /* x**3 */

#define QUAD(x) x * x * x * x   /* x**4 */



#define GMAX(A, B)      ((A) > (B) ? (A) : (B)) /* assign maximum of a and b */

#define GMIN(A, B)      ((A) < (B) ? (A) : (B)) /* assign minimum of a and b */



#define IMOD(A, B)      (A) - (((A) / (B)) * (B)) /* Integer mod function */



/* Projection codes



   0 = Geographic

   1 = Universal Transverse Mercator (UTM)

   2 = State Plane Coordinates

   3 = Albers Conical Equal Area

   4 = Lambert Conformal Conic

   5 = Mercator

   6 = Polar Stereographic

   7 = Polyconic

   8 = Equidistant Conic

   9 = Transverse Mercator

  10 = Stereographic

  11 = Lambert Azimuthal Equal Area

  12 = Azimuthal Equidistant

  13 = Gnomonic

  14 = Orthographic

  15 = General Vertical Near-Side Perspective

  16 = Sinusoidal

  17 = Equirectangular

  18 = Miller Cylindrical

  19 = Van der Grinten

  20 = (Hotine) Oblique Mercator 

  21 = Robinson

  22 = Space Oblique Mercator (SOM)

  23 = Alaska Conformal

  24 = Interrupted Goode Homolosine 

  25 = Mollweide

  26 = Interrupted Mollweide

  27 = Hammer

  28 = Wagner IV

  29 = Wagner VII

  30 = Oblated Equal Area

  99 = User defined

*/



/* Define projection codes */

enum ProjCode 

{

	NONE=-1,

	GEO=0, 

	_UTM, 

	SPCS, 

	ALBERS, 

	LAMCC, 

	MERCAT, 

	PS, 

	POLYC, 

	EQUIDC,

	TM, 

	STEREO, 

	LAMAZ, 

	AZMEQD,

	GNOMON, 

	ORTHO, 

	GVNSP, 

	SNSOID, 

	EQRECT, 

	MILLER, 

	VGRINT, 

	HOM, 

	ROBIN, 

	SOM, 

	ALASKA, 

	GOOD,

	MOLL, 

	IMOLL, 

	HAMMER, 

	WAGIV, 

	WAGVII, 

	OBEQA, 

	USDEF=99

};



/* Define unit code numbers to their names */

enum ProjUnit {

	  UNDEF=-1,

	  RADIAN=0, //Radians

	  FEET, //Feet

	  METER, //Meters

	  SECOND, //Seconds

	  DEGREE, //Decimal Degrees

	  INT_FEET //International Feet

};



//Spheroid codes

enum ProjDatum {

	NOT_SET=-1,

	CLARKE_1866 = 0,

	CLARKE_1880,

	BESSEL,

	INTERNATIONAL_1967,

	INTERNATIONAL_1909,

	WGS_72,

	EVEREST,

	WGS_66,

	GRS_1980,

	AIRY,

	MODIFIED_EVEREST,

	MODIFIED_AIRY,

	WGS_84,

	SOUTHEAST_ASIA,

	AUSTRALIAN_NATIONAL,

	KRASSOVSKY,

	HOUGH,

	MERCURY_1960,

	MODIFIED_MERCURY_1968,

	EARTH,

	BESSEL_1841_NAMIBIA,

	EVEREST_SABAH,

	EVEREST_INDIA_1956,

	EVEREST_MALAYSIA_1969,

	EVEREST_MALAY_1948,

	EVEREST_PAKISTAN,

	HAYFORD,

	HELMERT_1906,

	INDONESIAN_1974,

	SOUTH_AMERICAN_1969,

	WGS_60

};



//unit conversion factors

extern double factors[6][6];



extern long NADUT[134];



/* the Nad 27 State Plane Zones are set in this array

  --------------------------------------------------*/

extern long NAD27[134];



/* the Nad 83 State Plane Zones are set in this array

  --------------------------------------------------*/



extern long NAD83[134];



/* Semi-Major axis of supported Spheroids */

extern double major[SPHDCT];





/* Semi-Minor axis of supported Spheroids */

extern double minor[SPHDCT];

#endif

