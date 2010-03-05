
#include "hammer.h"

Hammer::Hammer() : Projection() 
{
	setNumber(HAMMER);
	setName("Hammer");
}

Hammer::Hammer(double gctpParams[], ProjUnit units, ProjDatum dat) :
Projection(gctpParams, units, dat) 
{
	setNumber(HAMMER);
	setName("Hammer");
}

void Hammer::_init() 
{
	return;
}

void Hammer::_inverse(double x, double y)
{

	double fac;

	/* Inverse equations
	-----------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;
	fac = sqrt(4.0 * m_radius * m_radius - (x * x)/ 4.0 - y * y) / 2.0;
	m_longitude = Util::adjust_lon(m_centerLon + 2.0 * 
					atan2((x * fac), (2.0 * m_radius * m_radius - x * x/4 - y * y)));
	m_latitude = Util::asinz(y * fac / m_radius / m_radius);
}

void Hammer::_forward(double lon, double lat)
{
	double dlon;
	double fac;

	/* Forward equations
	-----------------*/
	dlon = Util::adjust_lon(lon - m_centerLon);

	fac  = m_radius * 1.414213562 / sqrt(1.0 + cos(lat) * cos(dlon / 2.0));
	m_x_coord = m_falseEasting + fac * 2.0 * cos(lat) * sin(dlon / 2.0);
	m_y_coord = m_falseNorthing + fac * sin(lat);
}
