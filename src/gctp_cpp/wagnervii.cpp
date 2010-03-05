
#include "wagnervii.h"

WagnerVII::WagnerVII() : Projection()
{
	setNumber(WAGVII);
	setName("Wagner VII");
}

WagnerVII::WagnerVII(double gctpParams[], ProjUnit units, ProjDatum dat) :
Projection(gctpParams, units, dat)
{
	setNumber(WAGVII);
	setName("Wagner VII");
}

void WagnerVII::_init()
{
	return;
}

void WagnerVII::_forward(double lon, double lat)
{
	double delta_lon;	/* Delta longitude (Given longitude - center */
	double sin_lon, cos_lon;
	double s, c0, c1;

	/* Forward equations
	-----------------*/
	delta_lon = Util::adjust_lon(lon - m_centerLon);
	Util::gctp_sincos((delta_lon/3.0), &sin_lon, &cos_lon);
	s = 0.90631 * sin(lat);
	c0 = sqrt(1-s*s);
	c1 = sqrt(2.0 / (1.0 + c0 * cos_lon));
	m_x_coord = 2.66723 * m_radius * c0 * c1 * sin_lon + m_falseEasting;
	m_y_coord = 1.24104 * m_radius * s * c1 + m_falseNorthing;
}

void WagnerVII::_inverse(double x, double y)
{
	double t1, t2, p, c;

	/* Inverse equations
	-----------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;
	t1 = x / 2.66723;
	t2 = y / 1.24104;
	t1 *= t1;
	t2 *= t2;
	p = sqrt(t1 + t2);
	c = 2.0 * Util::asinz(p / (2.0 * m_radius));
	m_latitude = Util::asinz(y * sin(c) / (1.24104 * 0.90631 * p));
	m_longitude = Util::adjust_lon(m_centerLon + 3.0 * atan2(x * tan(c), 2.66723 * p));
}
