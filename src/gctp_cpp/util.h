#ifndef UTIL_H 
#define UTIL_H

#include "constants.h"
#include <math.h>


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
S. Nelson, EROS		June, 1993	Added inline commentswow,
S. Nelson, EROS		Nov, 1993	Added loop counter in ADJUST_LON
S. Nelson, EROS		Jan, 1998	Changed misspelled error message

*******************************************************************************/

class Util {
public:

	static void gctp_sincos(double val, double* sin_val, double* cos_val);
	static double asinz(double con);
	static double msfnz(double eccent, double sinphi, double cosphi);
	static double qsfnz(double eccent, double sinphi);
	static double phi1z(double eccent, double qs, long* flag);
	static double phi2z(double eccent, double ts, long* flag);
	static double phi3z(double ml, double e0, double e1, double e2, double e3, long* flag);
	static long phi4z(double eccent, double e0, double e1, double e2, double e3, double a, double b, double* c, double* phi);
	static double pakcz(double pak);
	static double paksz(double ang, long* iflg);
	static double pakr2dm(double pak);
	
	static double tsfnz(double eccent, double phi, double sinphi);
	static int sign(double x);
	static double adjust_lon(double x);

	static double e0fn(double x) {return(1.0-0.25*x*(1.0+x/16.0*(3.0+1.25*x)));}
	static double e1fn(double x) {return(0.375*x*(1.0+0.25*x*(1.0+0.46875*x)));}
	static double e2fn(double x) {return(0.05859375*x*x*(1.0+0.75*x));}
	static double e3fn(double x) {return(x*x*x*(35.0/3072.0));}
	static double e4fn(double x) { double con;
								   double com;
								   con = 1.0 + x;
								   com = 1.0 - x;
								   return (sqrt((pow(con,con))*(pow(com,com))));
								 }
	static double mlfn(double e0, double e1, double e2, double e3, double phi) {return(e0*phi-e1*sin(2.0*phi)+e2*sin(4.0*phi)-e3*sin(6.0*phi));}

	static long calc_utm_zone(double lon) {return((long)(((lon + 180.0) / 6.0) + 1.0));}
	
	static long untfz(long inunit, long outunit, double* factor);

	static void sphdz(long isph, double* parm, double* r_major, double* r_minor, double* radius);

	static long convertCoords(int fromUnit, int toUnit, double& x, double& y);

	//convert a packed DMS angle to radians
	static long DMSToRad(double& angle) 
	{
		long err = 0;
		angle = paksz(angle, &err) * 3600 * S2R;
		return(err);
	}
};

#endif




