
#include "orthographic.h"

Orthographic::Orthographic() : Projection(), m_sinCenterLat(0.0), m_cosCenterLat(0.0)
{
	setNumber(ORTHO);
	setName("Orthographic");
}

Orthographic::Orthographic(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat), m_sinCenterLat(0.0), m_cosCenterLat(0.0)
{
	setNumber(ORTHO);
	setName("Orthographic");
}

void Orthographic::_init() 
{
	Util::gctp_sincos(m_centerLat,&m_sinCenterLat,&m_cosCenterLat);
}

void Orthographic::_forward(double lon, double lat)
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
	ksp = 1.0;
	if ((g > 0) || (fabs(g) <= EPSLN))
	{
		m_x_coord = m_falseEasting + m_rMajor * ksp * cosphi * sin(dlon);
		m_y_coord = m_falseNorthing + m_rMajor * ksp * (m_cosCenterLat * sinphi -
					m_sinCenterLat * cosphi * coslon);
	}
	else
	{
		setError(143);
		return;
	}
}

void Orthographic::_inverse(double x, double y)
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
	if (rh > m_rMajor + .0000001)
	{
		setError(145);
		return;
	}
	z = Util::asinz(rh / m_rMajor);
	Util::gctp_sincos(z,&sinz,&cosz);
	m_longitude = m_centerLon;
	if (fabs(rh) <= EPSLN)
	{
		m_latitude = m_centerLat;
		return;
	}
	m_latitude = Util::asinz(cosz * m_sinCenterLat + (y * sinz * m_cosCenterLat)/rh);
	con = fabs(m_centerLat) - HALF_PI;
	if (fabs(con) <= EPSLN)
	{
		if (m_centerLat >= 0)
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
	con = cosz - m_sinCenterLat * sin(m_latitude);
	if ((fabs(con) >= EPSLN) || (fabs(x) >= EPSLN))
		m_longitude = Util::adjust_lon(m_centerLon + atan2((x * sinz * m_cosCenterLat), (con * rh)));
}