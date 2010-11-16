#ifndef CPROJ_H 
#define CPROJ_H


#include <math.h>
//#include "proj.h"

#define PI 	3.141592653589793238
#define HALF_PI (PI*0.5)
#define TWO_PI 	(PI*2.0)
#define EPSLN	1.0e-10
#define R2D     57.2957795131
/*
#define D2R     0.0174532925199
*/
#define D2R     1.745329251994328e-2
#define S2R	4.848136811095359e-6

#define OK	0
#define ERROR  -1
#define IN_BREAK -2


#define NULL 0
/* Misc macros
  -----------*/
//#define SQUARE(x)       x * x   /* x**2 */
//#define CUBE(x)     x * x * x   /* x**3 */
//#define QUAD(x) x * x * x * x   /* x**4 */

//#define GMAX(A, B)      ((A) > (B) ? (A) : (B)) /* assign maximum of a and b */
//#define GMIN(A, B)      ((A) < (B) ? (A) : (B)) /* assign minimum of a and b */

//#define IMOD(A, B)      (A) - (((A) / (B)) * (B)) /* Integer mod function */




#ifndef CPROJ_C
#define CPROJ_C

/*******************************************************************************
NAME                Projection support routines listed below.

PURPOSE:	The following functions are included in CPROJ.C.

		SINCOS:	  Calculates the sine and cosine.
		ASINZ:	  Eliminates roundoff errors.
		MSFNZ:	  Computes the constant small m for Oblique Equal Area.
		QSFNZ:	  Computes the constant small q for Oblique Equal Area.
		PHI1Z:	  Computes phi1 for Albers Conical Equal-Area.
		PHI2Z:	  Computes the latitude angle, phi2, for Lambert
			  Conformal Conic and Polar Stereographic.
		PHI3Z:	  Computes the latitude, phi3, for Equidistant Conic.
		PHI4Z:	  Computes the latitude, phi4, for Polyconic.
		PAKCZ:	  Converts a 2 digit alternate packed DMS format to
			  standard packed DMS format.
		PAKR2DM:  Converts radians to 3 digit packed DMS format.
		TSFNZ:	  Computes the small t for Lambert Conformal Conic and
			  Polar Stereographic.
		SIGN:	  Returns the sign of an argument.
		ADJUST_LON:  Adjusts a longitude angle to range -180 to 180.
		E0FN, E1FN, E2FN, E3FN:
			  Computes the constants e0,e1,e2,and e3 for
			  calculating the distance along a meridian.
		E4FN:	  Computes e4 used for Polar Stereographic.
		MLFN:	  Computes M, the distance along a meridian.
		CALC_UTM_ZONE:	Calculates the UTM zone number.

PROGRAMMER              DATE		REASON
----------              ----		------
D. Steinwand, EROS      July, 1991	Initial development
T. Mittan, EROS		May, 1993	Modified from Fortran GCTP for C GCTP
S. Nelson, EROS		June, 1993	Added inline comments
S. Nelson, EROS		Nov, 1993	Added loop counter in ADJUST_LON
S. Nelson, EROS		Jan, 1998	Changed misspelled error message

*******************************************************************************/

#define MAX_VAL 4
#define MAXLONG 2147483647.
#define DBLLONG 4.61168601e18

/* Function to calculate the sine and cosine in one call.  Some computer
   systems have implemented this function, resulting in a faster implementation
   than calling each function separately.  It is provided here for those
   computer systems which don`t implement this function
  ----------------------------------------------------*/
void sincos(double val,double *sin_val,double *cos_val );

/* Function to eliminate roundoff errors in asin
----------------------------------------------*/
double asinz ( double con );

/* Function to return the sign of an argument
  ------------------------------------------*/

int sign(double x );

/* Function to adjust a longitude angle to range from -180 to 180 radians
   added if statments 
  -----------------------------------------------------------------------*/
double adjust_lon( double x );		/* Angle in radians			*/


#endif //#ifndef CPROJ_C
#endif //def CPROJ_H


