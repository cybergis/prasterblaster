
#include "lambertazimuth.h"

LambertAzimuthal::LambertAzimuthal() : Projection(), m_sinCenterLat(0.0), m_cosCenterLat(0.0)
{
	setNumber(LAMAZ);
	setName("Lambert Azimuthal");
}

LambertAzimuthal::LambertAzimuthal(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat), m_sinCenterLat(0.0), m_cosCenterLat(0.0)
{
	setNumber(LAMAZ);
	setName("Lambert Azimuthal");
}

void LambertAzimuthal::_init()
{
	Util::gctp_sincos(m_centerLat, &m_sinCenterLat, &m_cosCenterLat);
}

void LambertAzimuthal::_forward(double lon, double lat)
{
	double delta_lon;	/* Delta longitude (Given longitude - center 	*/
	double sin_delta_lon;	/* Sine of the delta longitude 			*/
	double cos_delta_lon;	/* Cosine of the delta longitude 		*/
	double sin_lat;		/* Sine of the given latitude 			*/
	double cos_lat;		/* Cosine of the given latitude 		*/
	double g;		/* temporary varialbe				*/
	double ksp;		/* heigth above elipsiod			*/

	/* Forward equations
	-----------------*/
	delta_lon = Util::adjust_lon(lon - m_centerLon);
	Util::gctp_sincos(lat, &sin_lat, &cos_lat);
	Util::gctp_sincos(delta_lon, &sin_delta_lon, &cos_delta_lon);
	g = m_sinCenterLat * sin_lat + m_cosCenterLat * cos_lat * cos_delta_lon;
	if (g == -1.0) 
	{
		setError(113);
		return;
	}

	ksp = m_radius * sqrt(2.0 / (1.0 + g));
	m_x_coord = ksp * cos_lat * sin_delta_lon + m_falseEasting;
	m_y_coord = ksp * (m_cosCenterLat * sin_lat - m_sinCenterLat * cos_lat * m_cosCenterLat) + 
				m_falseNorthing;
}

void LambertAzimuthal::_inverse(double x, double y)
{
	double Rh;
	double z;		/* Great circle dist from proj center to given point */
	double sin_z;		/* Sine of z */
	double cos_z;		/* Cosine of z */
	double temp;		/* Re-used temporary variable */


	/* Inverse equations
	-----------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;
	Rh = sqrt(x * x + y * y);
	temp = Rh / (2.0 * m_radius);
	if (temp > 1) 
	{
		setError(115);
		return;
	}

	z = 2.0 * Util::asinz(temp);
	Util::gctp_sincos(z, &sin_z, &cos_z);
	m_longitude = m_centerLon;
	if (fabs(Rh) > EPSLN)
	{
		m_latitude = Util::asinz(m_sinCenterLat * cos_z + m_cosCenterLat * sin_z * y / Rh);
		temp = fabs(m_centerLat) - HALF_PI;
		if (fabs(temp) > EPSLN)
		{
			temp = cos_z - m_sinCenterLat * sin(m_latitude);

			if(temp!=0.0)
				m_longitude = Util::adjust_lon(m_centerLon + atan2(x*sin_z*m_cosCenterLat,temp*Rh));
		}
		else if (m_centerLat < 0.0) 
			m_longitude = Util::adjust_lon(m_centerLon - atan2(-x, y));
		else 
			m_longitude = Util::adjust_lon(m_centerLon + atan2(x, -y));
	}
	else 
		m_latitude = m_centerLat;
}