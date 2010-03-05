
#include "miller.h"

Miller::Miller() : Projection()
{
	setNumber(MILLER);
	setName("Miller_Cylindrical");
}

Miller::Miller(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat)
{
	setNumber(MILLER);
	setName("Miller_Cylindrical");
}

void Miller::_init() 
{
	return;
}

void Miller::_forward(double lon, double lat)
{
	double dlon;

	/* Forward equations
	-----------------*/
	dlon = Util::adjust_lon(lon - m_centerLon);
	m_x_coord = m_falseEasting + m_radius * dlon;
	m_y_coord = m_falseNorthing + m_radius * log(tan((PI / 4.0) + (lat / 2.5))) * 1.25;

}

void Miller::_inverse(double x, double y)
{

	/* Inverse  equations
	------------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;

	m_longitude = Util::adjust_lon(m_centerLon + x / m_radius);
	m_latitude = 2.5 * (atan(exp(y / m_radius / 1.25)) - PI / 4.0);


}

