#include "cproj.h"

void sincos(double val,double *sin_val,double *cos_val )
{
 *sin_val = sin(val);
 *cos_val = cos(val);
 return;
}

double asinz ( double con )
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

double adjust_lon( double x )		/* Angle in radians			*/
{
long count = 0;
for(;;)
  {
  if (fabs(x)<=PI)
     break;
  else
  if (((long) fabs(x / PI)) < 2)
     x = x-(sign(x) *TWO_PI);
  else
  if (((long) fabs(x / TWO_PI)) < MAXLONG)
     {
     x = x-(((long)(x / TWO_PI))*TWO_PI);
     }
  else
  if (((long) fabs(x / (MAXLONG * TWO_PI))) < MAXLONG)
     {
     x = x-(((long)(x / (MAXLONG * TWO_PI))) * (TWO_PI * MAXLONG));
     }
  else
  if (((long) fabs(x / (DBLLONG * TWO_PI))) < MAXLONG)
     {
     x = x-(((long)(x / (DBLLONG * TWO_PI))) * (TWO_PI * DBLLONG));
     }
  else
     x = x-(sign(x) *TWO_PI);
  count++;
  if (count > MAX_VAL)
     break;
  }

return(x);
}



int sign(double x )
{
if (x < 0.0)
    return(-1);
else
    return(1);
}

