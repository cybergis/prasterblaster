
#include "gnomonic.h"

Gnomonic::Gnomonic() : Projection()
{
	setNumber(GNOMON);
	setName("Gnomonic");
}

Gnomonic::Gnomonic(double gctpParams[], ProjUnit units, ProjDatum dat) :
Projection(gctpParams, units, dat)
{
	setNumber(GNOMON);
	setName("Gnomonic");
}

void Gnomonic::_init()
{
	Util::gctp_sincos(m_centerLat, &m_sinCenterLat, &m_cosCenterLat);
}

void Gnomonic::_forward(double lon, double lat)
{
	double dlon;
	double sinphi,cosphi;
	double coslon;
	double g;
	double ksp;


	/* Forward equations
	-----------------*/
	dlon = Util::adjust_lon(lon - m_centerLon);
	Util::gctp_sincos(lat,&sinphi,&cosphi);
	coslon = cos(dlon);
	g = m_sinCenterLat * sinphi + m_cosCenterLat * cosphi * coslon;
	if (g <= 0.0)
	{
		setError(133);
		return;
	}
	ksp = 1.0 / g;
	m_x_coord = m_falseEasting + m_radius * ksp * cosphi * sin(dlon);
	m_y_coord = m_falseNorthing + m_radius * ksp * (m_cosCenterLat * sinphi - m_sinCenterLat * cosphi * 
				coslon);
}

void Gnomonic::_inverse(double x, double y)
{

	double rh;
	double z,sinz,cosz;
	double con;


	/* Inverse equations
	-----------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;
	rh = sqrt(x * x + y * y);
	z = atan(rh / m_radius);
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
	if ((fabs(con) < EPSLN) && (fabs(x) < EPSLN))
		return;

	m_longitude = Util::adjust_lon(m_centerLon + atan2((x * sinz * m_cosCenterLat), (con * rh)));
}