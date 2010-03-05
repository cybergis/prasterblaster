
#include "stereo.h"

Stereo::Stereo() : Projection(), m_sinCenterLat(0.0), m_cosCenterLat(0.0)
{
	setNumber(STEREO);
	setName("Stereographic");
}

Stereo::Stereo(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat), m_sinCenterLat(0.0), m_cosCenterLat(0.0)
{
	setNumber(STEREO);
	setName("Stereographic");
}

void Stereo::_init()
{
	Util::gctp_sincos(m_centerLat,&m_sinCenterLat,&m_cosCenterLat);
}

void Stereo::_forward(double lon, double lat)
{

	double sinphi, cosphi;	/* sin and cos value				*/
	double dlon;		/* delta longitude value			*/
	double coslon;		/* cos of longitude				*/
	double ksp;		/* scale factor					*/
	double g;		

	/* Forward equations
	-----------------*/
	dlon = Util::adjust_lon(lon - m_centerLon);
	Util::gctp_sincos(lat,&sinphi,&cosphi);
	coslon = cos(dlon);
	g = m_sinCenterLat * sinphi + m_cosCenterLat * cosphi * coslon;

	if (fabs(g + 1.0) <= EPSLN)
	{
		setError(103);
		return;
	}

	ksp = 2.0 / (1.0 + g);
	m_x_coord = m_falseEasting + m_rMajor * ksp * cosphi * sin(dlon);
	m_y_coord = m_falseNorthing + m_rMajor * ksp * (m_cosCenterLat * sinphi - m_sinCenterLat * 
				cosphi * coslon);
}

void Stereo::_inverse(double x, double y)
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
	z = 2.0 * atan(rh / (2.0 * m_rMajor));
	
	Util::gctp_sincos(z,&sinz,&cosz);
	
	m_longitude = m_centerLon;
	if (fabs(rh) <= EPSLN)
	{
		m_latitude = m_centerLat;
		return;
	}
	else
	{
		m_latitude = asin(cosz * m_sinCenterLat + (y * sinz * m_cosCenterLat) / rh);
		con = fabs(m_centerLat) - HALF_PI;
		if (fabs(con) <= EPSLN)
		{
			if (m_centerLat >= 0.0)
			{
				m_longitude = Util::adjust_lon(m_centerLon + atan2(x, -y));
				return;
			}
			else
			{
				m_longitude = Util::adjust_lon(m_centerLon - atan2(-x, y));
				return;
			}
		}
		else
		{
			con = cosz - m_sinCenterLat * sin(m_latitude);
			if ((fabs(con) < EPSLN) && (fabs(x) < EPSLN))
				return;
			else
				m_longitude = Util::adjust_lon(m_centerLon + 
											   atan2((x * sinz * m_cosCenterLat), (con * rh)));
		}
	}
}
	