
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
#define NULL 0
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
static double factors[6][6] = {
	{1.0, 0.0, 0.0, 206264.8062470963, 57.29577951308231, 0.0},
	{0.0, 1.0, .3048006096012192, 0.0, 0.0, 1.000002000004},
	{0.0, 3.280833333333333, 1.0, 0.0, 0.0, 3.280839895013124},
	{.484813681109536e-5, 0.0, 0.0, 1.0, .27777777777778e-3, 0.0}, 
	{.0174532925199433, 0.0, 0.0, 3600, 1.0, 0.0},
	{0.0, .999998, .3048, 0.0, 0.0, 1.0}
};

static long NADUT[134] = {1, 5, 1, 1, 5, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2,
			  1, 1, 5, 2, 1, 2, 5, 1, 2, 2, 2, 1, 1, 1, 5, 2, 1, 5,
			  2, 2, 5, 2, 1, 1, 5, 2, 2, 1, 2, 1, 2, 2, 1, 2, 2, 2
};
/* the Nad 27 State Plane Zones are set in this array
  --------------------------------------------------*/
static long NAD27[134] = {101,102,5010,5300,201,202,203,301,302,401,402,403,404,
				405,406,407,501,502,503,600,700,901,902,903,1001,1002,5101,
				5102,5103,5104,5105,1101,1102,1103,1201,1202,1301,1302,1401,
				1402,1501,1502,1601,1602,1701,1702,1703,1801,1802,1900,2001,
				2002,2101,2102,2103,2111,2112,2113,2201,2202,2203,2301,2302,
				2401,2402,2403,2501,2502,2503,2601,2602,2701,2702,2703,2800,
				2900,3001,3002,3003,3101,3102,3103,3104,3200,3301,3302,3401,
				3402,3501,3502,3601,3602,3701,3702,3800,3901,3902,4001,4002, 
				4100,4201,4202,4203,4204,4205,4301,4302,4303,4400,4501,4502,
				4601,4602,4701,4702,4801,4802,4803,4901,4902,4903,4904,5001,
				5002,5003,5004,5005,5006,5007,5008,5009,5201,5202,5400};

/* the Nad 83 State Plane Zones are set in this array
  --------------------------------------------------*/

static long NAD83[134] = {101,102,5010,5300,201,202,203,301,302,401,402,403,
                404,405,406,0000,501,502,503,600,700,901,902,903,1001,1002,
                5101,5102,5103,5104,5105,1101,1102,1103,1201,1202,1301,1302,
                1401,1402,1501,1502,1601,1602,1701,1702,1703,1801,1802,1900,
                2001,2002,2101,2102,2103,2111,2112,2113,2201,2202,2203,2301,
                2302,2401,2402,2403,2500,0000,0000,2600,0000,2701,2702,2703,
                2800,2900,3001,3002,3003,3101,3102,3103,3104,3200,3301,3302,
                3401,3402,3501,3502,3601,3602,3701,3702,3800,3900,0000,4001,
                4002,4100,4201,4202,4203,4204,4205,4301,4302,4303,4400,4501,
                4502,4601,4602,4701,4702,4801,4802,4803,4901,4902,4903,4904,
                5001,5002,5003,5004,5005,5006,5007,5008,5009,5200,0000,5400
};

	/* Semi-Major axis of supported Spheroids */
static double major[SPHDCT] = {
		6378206.4,		/* 0: Clarke 1866  */
		6378249.145,		/* 1: Clarke 1880 */
		6377397.155,		/* 2: Bessel */
		6378157.5,		/* 3: International 1967 */
		6378388.0,		/* 4: International 1909 */
		6378135.0,		/* 5: WGS 72 */
		6377276.3452,		/* 6: Everest */
		6378145.0,		/* 7: WGS 66 */
        6378137.0,		/* 8: GRS 1980 */
		6377563.396,		/* 9: Airy */
		6377304.063,		/* 10: Modified Everest */
		6377340.189,		/* 11: Modified Airy */
        6378137.0,		/* 12: WGS 84 */
		6378155.0,		/* 13: Southeast Asia */
		6378160.0,		/* 14: Australian National */
		6378245.0,		/* 15: Krassovsky */
        6378270.0,		/* 16: Hough */
		6378166.0,		/* 17: Mercury 1960 */
		6378150.0,		/* 18: Modified Mercury 1968 */
		6370997.0,		/* 19: Sphere of Radius 6370997 meters*/
		6377483.865,		/* 20: Bessel 1841(Namibia) */
		6377298.556,		/* 21: Everest (Sabah & Sarawak) */
		6377301.243,		/* 22: Everest (India 1956) */
		6377295.664,		/* 23: Everest (Malaysia 1969) */
		6377304.063,		/* 24: Everest (Malay & Singapr 1948)*/
		6377309.613,		/* 25: Everest (Pakistan) */
		6378388.0,		/* 26: Hayford */
		6378200.0,		/* 27: Helmert 1906 */
		6378160.000,		/* 28: Indonesian 1974 */
		6378160.0,		/* 29: South American 1969 */
		6378165.0
};		/* 30: WGS 60 */

	/* Semi-Minor axis of supported Spheroids */
static double minor[SPHDCT] = {
		6356583.8,		/* 0: Clarke 1866 */
		6356514.86955,		/* 1: Clarke 1880 */
		6356078.96284,		/* 2: Bessel */
		6356772.2,		/* 3: International 1967 */
        6356911.94613,		/* 4: International 1909 */
		6356750.519915,		/* 5: WGS 72 */
		6356075.4133,		/* 6: Everest */
        6356759.769356,		/* 7: WGS 66 */
		6356752.31414,		/* 8: GRS 1980 */
		6356256.91,		/* 9: Airy */
        6356103.039,		/* 10: Modified Everest */
		6356034.448,		/* 11: Modified Airy */
		6356752.314245,		/* 12: WGS 84 */
        6356773.3205,		/* 13: Southeast Asia */
		6356774.719,		/* 14: Australian National */
		6356863.0188,		/* 15: Krassovsky */
        6356794.343479,		/* 16: Hough */
		6356784.283666,		/* 17: Mercury 1960 */
		6356768.337303,		/* 18: Modified Mercury 1968 */
        6370997.0,		/* 19: Sphere of Radius 6370997 meters*/
		6356165.382966,		/* 20: Bessel 1841(Namibia) */
		6356097.571445,		/* 21: Everest (Sabah & Sarawak) */
		6356100.228368,		/* 22: Everest (India 1956) */
		6356094.667915,		/* 23: Everest (Malaysia 1969) */
		6356103.038993,		/* 24: Everest (Malay & Singapr 1948)*/
		6356108.570542,		/* 25: Everest (Pakistan) */
		6356911.946128,		/* 26: Hayford */
		6356818.169,		/* 27: Helmert 1906 */
		6356774.504086,		/* 28: Indonesian 1974 */
		6356774.719,		/* 29: South American 1969 */
		6356783.287
};		/* 30: WGS 60 */

#endif