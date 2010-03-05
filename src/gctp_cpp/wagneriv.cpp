#include "wagneriv.h"


WagnerIV::WagnerIV() : Projection()
{
	setNumber(WAGIV);
	setName("Wagner IV");
}

WagnerIV::WagnerIV(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat)
{
	setNumber(WAGIV);
	setName("Wagner IV");
}

void WagnerIV::_init()
{
	return;
}

void WagnerIV::_forward(double lon, double lat)
{

	double delta_lon;	/* Delta longitude (Given longitude - center */
	double theta;
	double delta_theta;
	double con;
	long i;

	/* Forward equations
	-----------------*/
	delta_lon = Util::adjust_lon(lon - m_centerLon);
	theta = lat;
	con = 2.9604205062 * sin(lat);

	/* Iterate using the Newton-Raphson method to find theta
	-----------------------------------------------------*/
	for (i=0;;i++)
	{
		delta_theta = -(theta + sin(theta) - con) / (1.0 + cos(theta));
		theta += delta_theta;
		if (fabs(delta_theta) < EPSLN) 
			break;
		if (i >= 30) {
			setError(2);
			return;
		}
	}
	theta /= 2.0;
	m_x_coord = 0.86310 * m_radius * delta_lon * cos(theta) + m_falseEasting;
	m_y_coord = 1.56548 * m_radius * sin(theta) + m_falseNorthing;
}

void WagnerIV::_inverse(double x, double y)
{
	double theta;

	/* Inverse equations
	-----------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;
	theta = asin(y /  (1.56548 * m_radius));
	m_longitude = Util::adjust_lon(m_centerLon + (x / (0.86310 * m_radius * cos(theta))));
	m_latitude = asin((2.0 * theta + sin(2.0 * theta)) / 2.9604205062);
}