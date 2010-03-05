
#include "azequidistant.h"

AzEquidistant::AzEquidistant() : Projection(), m_sinCenterLat(0.0), m_cosCenterLat(0.0)
{
   setNumber(AZMEQD);
   setName("Azimuthal Equidistant");
}

AzEquidistant::AzEquidistant(double gctpParams[], ProjUnit units, ProjDatum dat) :
Projection(gctpParams, units, dat), m_sinCenterLat(0.0), m_cosCenterLat(0.0)
{
   setNumber(AZMEQD);
   setName("Azimuthal Equidistant");
}

void AzEquidistant::_init()
{
   Util::gctp_sincos(m_centerLat,&m_sinCenterLat,&m_cosCenterLat);
}

void AzEquidistant::_forward(double lon, double lat)
{
   double sinphi, cosphi;	/* sin and cos value				*/
   double dlon;		/* delta longitude value			*/
   double coslon;		/* cos of longitude				*/
   double ksp;		/* scale factor					*/
   double g;		
   double z;		/* angle					*/

   dlon = Util::adjust_lon(lon - m_centerLon);
   Util::gctp_sincos(lat,&sinphi,&cosphi);
   coslon = cos(dlon);
   g = m_sinCenterLat * sinphi + m_cosCenterLat * cosphi * coslon;
   if (fabs(fabs(g) - 1.0) < EPSLN)
   {
      ksp = 1.0;
      if (g < 0.0)
      {
         //con = 2.0 * HALF_PI * r_major;
         //sprintf(mess,"Point projects into a circle of radius = %12.2lf",con);
         //p_error(mess,"azim-for");  
         setError(123);
         return;
      }
   }
   else
   {
      z = acos(g);
      ksp = z/ sin(z);
   }

   m_x_coord = m_falseEasting + m_rMajor * ksp * cosphi * sin(dlon);
   m_y_coord = m_falseNorthing + m_rMajor * ksp * (m_cosCenterLat * sinphi - m_sinCenterLat * 
      cosphi * coslon);

}

void AzEquidistant::_inverse(double x, double y)
{
   double rh;		/* height above ellipsoid			*/
   double z;		/* angle					*/
   double sinz,cosz;	/* sin of z and cos of z			*/
   double con;


   /* Inverse equations
   -----------------*/
   x -= m_falseEasting;
   y -= m_falseNorthing;
   rh = sqrt(x * x + y * y);
   if (rh > (2.0 * HALF_PI * m_rMajor))
   {
      setError(125);
      return;
   }

   z = rh / m_rMajor;
   Util::gctp_sincos(z,&sinz,&cosz);
   m_longitude = m_centerLon;
   if (fabs(rh) <= EPSLN)
   {
      m_latitude = m_centerLat;
      return;
   }

   m_latitude = Util::asinz(cosz * m_sinCenterLat + (y * sinz * m_cosCenterLat) / rh);
   con = fabs(m_centerLat) - HALF_PI;
   if (fabs(con) <= EPSLN)
   {
      if (m_centerLat >= 0.0)
      {
         m_longitude = Util::adjust_lon(m_centerLon + atan2(x , -y));
         return;
      }
      else
      {
         m_longitude = Util::adjust_lon(m_centerLon - atan2(-x , y));
         return;
      }
   }

   con = cosz - m_sinCenterLat * sin(m_latitude);
   if ((fabs(con) < EPSLN) && (fabs(x) < EPSLN))
   {
      return;
   }

  
   m_longitude = Util::adjust_lon(m_centerLon + atan2((x * sinz * m_cosCenterLat), (con * rh)));

}