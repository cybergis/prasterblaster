#include "util.h"
#include <stdio.h>


void Util::gctp_sincos(double val,double *sin_val,double *cos_val )
{
 *sin_val = sin(val);
 *cos_val = cos(val);
 return;
}

double Util::asinz ( double con )
{
 if (fabs(con) > 1.0)
   {
   if (con > 1.0)
     con = 1.0;
   else
     con = -1.0;
   }
 return(asin(con));
}

double Util::adjust_lon( double x )		/* Angle in radians			*/
{
	long count = 0;
	for(;;)
	{
		
		if (fabs(x)<=PI)
			break;
		
		else if (((long) fabs(x / PI)) < 2)
			x = x-(sign(x) *TWO_PI);
		
		else if (((long) fabs(x / TWO_PI)) < MAXLONG)
			x = x-(((long)(x / TWO_PI))*TWO_PI);
		
		
		else if (((long) fabs(x / (MAXLONG * TWO_PI))) < MAXLONG)
			x = x-(((long)(x / (MAXLONG * TWO_PI))) * (TWO_PI * MAXLONG));
		
		
		else if (((long) fabs(x / (DBLLONG * TWO_PI))) < MAXLONG)
			x = x-(((long)(x / (DBLLONG * TWO_PI))) * (TWO_PI * DBLLONG));
		
		else
			x = x-(sign(x) *TWO_PI);
		
		count++;
		
		if (count > MAX_VAL)
			break;
	}

	return(x);
}


int Util::sign(double x )
{
if (x < 0.0)
    return(-1);
else
    return(1);
}

double Util::msfnz ( double eccent, double sinphi,   double cosphi )
{
	double con;

	con = eccent * sinphi;
	return((cosphi / (sqrt (1.0 - con * con))));
}

double Util::qsfnz (double eccent, double sinphi)
{
	double con;

	if (eccent > 1.0e-7)
	{
		con = eccent * sinphi;
		return (( 1.0- eccent * eccent) * (sinphi /(1.0 - con * con) - (.5/eccent)*
			log((1.0 - con)/(1.0 + con))));
	}
	else
		return(2.0 * sinphi);
}

double Util::phi1z ( double eccent,	/* Eccentricity angle in radians		*/
					  double qs,		/* Angle in radians				*/
			          long  *flag )	/* Error flag number				*/
{

	double eccnts;
	double dphi;
	double con;
	double com;
	double sinpi;
	double cospi;
	double phi;
	long i;

	phi = asinz(.5 * qs);

	if (eccent < EPSLN) 
		return(phi);
	
	eccnts = eccent * eccent;

	for (i = 1; i <= 25; i++) 
	{
		Util::gctp_sincos(phi,&sinpi,&cospi);
		con = eccent * sinpi; 
		com = 1.0 - con * con;
		dphi = .5 * com * com / cospi * (qs / (1.0 - eccnts) - sinpi / com + 
			   .5 / eccent * log ((1.0 - con) / (1.0 + con)));
		phi = phi + dphi;
		
		if (fabs(dphi) <= 1e-7)
			return(phi);
	}
	// p_error ("Convergence error","phi1z-conv");
	*flag = 001;
	return(ERROR);
}

double Util::phi2z(double eccent, double ts, long* flag) {
	double eccnth;
	double phi;
	double con;
	double dphi;
	double sinpi;
	long i;

	*flag = 0;
	eccnth = .5 * eccent;
	phi = HALF_PI - 2 * atan(ts);
	for (i = 0; i <= 15; i++)
	{
		sinpi = sin(phi);
		con = eccent * sinpi;
		dphi = HALF_PI - 2 * atan(ts *(pow(((1.0 - con)/(1.0 + con)),eccnth))) - 
			phi;
		phi += dphi; 
		if (fabs(dphi) <= .0000000001)
			return(phi);
	}
	*flag = 002;
	return(002);
}

double Util::phi3z(double ml, double e0, double e1, double e2, double e3, long* flag) {

	double phi;
	double dphi;
	long i;

	phi = ml;
	for (i = 0; i < 15; i++)
	{
		dphi = (ml + e1 * sin(2.0 * phi) - e2 * sin(4.0 * phi) + e3 * sin(6.0 * phi))
			/ e0 - phi;
		phi += dphi;
		if (fabs(dphi) <= .0000000001)
		{
			*flag = 0;
			return(phi);
		}
	}
	printf("Latitude failed to converge after 15 iterations -- PHI3Z-CONV");
	*flag = 3;
	return(3);
}

long Util::phi4z(double eccent, double e0, double e1, double e2, double e3, double a, double b, double* c, double* phi) {

	double sinphi;
	double sin2ph;
	double tanphi;
	double ml;
	double mlp;
	double con1;
	double con2;
	double con3;
	double dphi;
	long i;

	*phi = a;
	for (i = 1; i <= 15; i++)
	{
		sinphi = sin(*phi);
		tanphi = tan(*phi);
		*c = tanphi * sqrt (1.0 - eccent * sinphi * sinphi);
		sin2ph = sin (2.0 * *phi);
		ml = e0 * *phi - e1 * sin2ph + e2 * sin (4.0 *  *phi) - e3 * 
			sin (6.0 *  *phi);
		mlp = e0 - 2.0 * e1 * cos (2.0 *  *phi) + 4.0 * e2 *
			cos (4.0 *  *phi) - 6.0 * e3 * cos (6.0 *  *phi);
		con1 = 2.0 * ml + *c * (ml * ml + b) - 2.0 * a *  (*c * ml + 1.0);
		con2 = eccent * sin2ph * (ml * ml + b - 2.0 * a * ml) / (2.0 * *c);
		con3 = 2.0 * (a - ml) * (*c * mlp - 2.0 / sin2ph) - 2.0 * mlp;
		dphi = con1 / (con2 + con3);
		*phi += dphi;
		if (fabs(dphi) <= .0000000001 )
			return(OK);
	}
	printf("Latitude failed to converge -- phi4z-conv");
	return(004);
}

double Util::pakcz(double pak) {
	double con;
	double secs;
	long degs,mins;

	con = fabs (pak);
	degs = (long) ((con / 10000.0) + .001);
	con =  con  - degs * 10000;
	mins = (long) ((con / 100.0) + .001);
	secs = con  - mins * 100;
	con = (double) (degs) * 1000000.0 + (double) (mins) * 1000.0 + secs;

	if (pak < 0) 
		con *= -1;

	return(con); 
}

double Util::pakr2dm(double pak) {
  double con;
  double secs;
  long degs,mins;
  char sgna;

  sgna = ' ';
  pak *= R2D;

  if (pak < 0.0) 
	 sgna = '-';

  con = fabs (pak);
  degs = (long) (con);
  con =  (con  - degs) * 60;
  mins = (long) con;
  secs = (con  - mins) * 60;
  con = (double) (degs) * 1000000.0 + (double) (mins) * 1000.0 + secs;

  if (sgna == '-') 
	  con = - con;

  return(con); 
}

double Util::tsfnz(double eccent, double phi, double sinphi) {
  double con;
  double com;
  
  con = eccent * sinphi;
  com = .5 * eccent; 
  con = pow(((1.0 - con) / (1.0 + con)),com);
  return (tan(.5 * (HALF_PI - phi))/con);
}

long Util::untfz(long inunit, long outunit, double* factor) {
	if ((outunit >= 0) && (outunit <= MAXUNIT) && (inunit >= 0)
		&& (inunit <= MAXUNIT))
	{
		*factor = factors[inunit][outunit];

		/* Angle units can not be converted to length units
		------------------------------------------------*/
		if (*factor == 0.0)
		{
			printf("Incompatable unit codes -- untfz-code");
			return(1101);
		}
	}
	else
	{
		printf("Illegal source or target unit code -- untfz-unit\n");
		return(5);
	}

	return(OK);
}

double Util::paksz(double ang, long* iflg) {
	double fac;		/* sign flag			*/
	double deg;		/* degree variable		*/
	double min;		/* minute variable		*/
	double sec;		/* seconds variable		*/
	double tmp;		/* temporary variable		*/
	long i;			/* temporary variable		*/


	*iflg = 0;

	if (ang < 0.0)
		fac = -1;
	
	else
		fac = 1;

	/* find degrees
	-------------*/
	sec = fabs(ang);
	tmp = 1000000.0;
	i = (long) (sec/tmp);
	
	if (i > 360)
	{
		printf("Illegal DMS field -- paksz-deg");
		*iflg = 1116;
		return(ERROR);
	}
	
	else
		deg = i;

	/* find minutes
	-------------*/
	sec = sec - deg * tmp;
	tmp = 1000;
	i = (long) (sec/tmp);
	
	if (i > 60)
	{
		printf("Illegal DMS field -- paksz-min");
		*iflg = 1116;
		return(ERROR);
	}
	
	else
		min = i;

	/* find seconds
	-------------*/
	sec = sec - min * tmp;
	
	if (sec > 60)
	{
		printf("Illegal DMS field -- paksz-sec");
		*iflg = 1116;
		return(ERROR);
	}
	
	else
		sec = fac * (deg * 3600.0 + min * 60.0 + sec);
	
	deg = sec / 3600.0;

	return(deg);
}

void Util::sphdz(long isph, double* parm, double* r_major, double* r_minor, double* radius) {
	double t_major;		/* temporary major axis				*/
	double t_minor;		/* temporary minor axis				*/
	long jsph;		/* spheroid code number				*/

	/* if the spheroid code is a negative number, get the semi-major and semi-minor
	axis from the projection array
	--------------------------------------------------------------------------*/
	if (isph < 0)
	{
		t_major = fabs(parm[0]);
		t_minor = fabs(parm[1]);

		if (t_major  > 0.0) 
		{
			/* The semimajor axis and the semiminor axis are in the array, assign
			them directly
			--------------------------------------------------------------------*/
			if (t_minor > 1.0)
			{
				*r_major = t_major;
				*r_minor = t_minor;
				*radius = t_major;
			} 
			/* The semimajor axis and the eccentricity squared values are in the array,
			therefore, the semiminor axis is computed from the eccentricity
			squared value parm[1]
			----------------------------------------------------------------------*/
			else
				if (t_minor > 0.0)
				{
					*r_major = t_major;
					*radius = t_major;
					*r_minor = (sqrt(1.0 - t_minor)) * t_major; 
				}
				/* The semiminor axis is zero or less, assign the semimajor axis to
				the semiminor axis.
				-----------------------------------------------------------------*/
				else
				{
					*r_major = t_major;
					*radius = t_major;
					*r_minor = t_major;
				}
		}
		/* The sphroid code is to be used to assign the axis
		-------------------------------------------------*/
		else
			if (t_minor > 0.0)	/* t_major = 0 */

				/* Assign Clarke 1866 semi-major and semi-minor axis
				---------------------------------------------------*/
			{
				*r_major = major[0];
				*radius = major[0];
				*r_minor = minor[0];
			}
			else
				/* Assign Spheroid radius to semi-major and semi-minor axis
				---------------------------------------------------------*/
			{
				*r_major = major[RADVAL];
				*radius = major[RADVAL];
				*r_minor = major[RADVAL];
			}
	}
	/* Use the spheroid code to assign the semi-major and semi-minor axis
	-----------------------------------------------------------------*/
	else		/* isph >= 0 */
	{
		jsph = isph;

		/* The spheroid code is out of range, assign Clarke 1866
		------------------------------------------------------*/
		if (jsph > (SPHDCT - 1))
		{
			printf("Invalid spheroid selection -- INFORMATIONAL");
			printf("Reset to 0 -- INFORMATIONAL");
			jsph = 0;
		}
		/* Assign the radius argument to the standard radius value
		-------------------------------------------------------*/
		*r_major = major[jsph];
		*r_minor = minor[jsph];
		*radius = major[RADVAL];
	}
	return;
}

long Util::convertCoords(int fromUnit, int toUnit, double& x, double& y) {
	double factor = 0.0;
	long err = 0;
	err = Util::untfz(fromUnit, toUnit, &factor);
	
	x *= factor;
	y *= factor;
	return(err);
}









