#include "vandergrinten.h"

VanDerGrinten::VanDerGrinten() : Projection()
{
	setNumber(VGRINT);
	setName("Van Der Grinten");
}

VanDerGrinten::VanDerGrinten(double gctpParams[], ProjUnit units, ProjDatum dat) :
Projection(gctpParams, units, dat)
{
	setNumber(VGRINT);
	setName("Van Der Grinten");
	setParamLoad();
}


void VanDerGrinten::_init()
{
	return;
}

void VanDerGrinten::_forward(double lon, double lat)
{
	double dlon;
	double theta;
	double al,asq;
	double g,gsq;
	double m,msq;
	double con;
	double costh,sinth;

	/* Forward equations
	-----------------*/
	dlon = Util::adjust_lon(lon  - m_centerLon);

	if (fabs(lat) <= EPSLN)
	{
		m_x_coord = m_falseEasting  + m_radius * dlon;
		m_y_coord = m_falseNorthing;
		return;
	}

	theta = Util::asinz(2.0 * fabs(lat / PI));

	if ((fabs(dlon) <= EPSLN) || (fabs(fabs(lat) - HALF_PI) <= EPSLN))
	{
		m_x_coord = m_falseEasting;
		
		if (lat >= 0)
			m_y_coord = m_falseNorthing + PI * m_radius * tan(.5 * theta);
		else
			m_x_coord = m_falseNorthing + PI * m_radius * -tan(.5 * theta);

		return;
	}
	
	al = .5 * fabs((PI / dlon) - (dlon / PI));
	asq = al * al;
	Util::gctp_sincos(theta,&sinth,&costh);
	g = costh / (sinth + costh - 1.0);
	gsq = g * g;
	m = g * (2.0 / sinth - 1.0);
	msq = m * m;
	con = PI * m_radius * (al * (g - msq) + sqrt(asq * (g - msq) * (g - msq) - (msq + asq)
		* (gsq - msq))) / (msq + asq);
	
	if (dlon < 0)
		con = -con;
	
	m_x_coord = m_falseEasting + con;
	con = fabs(con / (PI * m_radius));

	if (lat >= 0)
		m_y_coord = m_falseNorthing + PI * m_radius * sqrt(1.0 - con * con - 2.0 * al * con);
	else
		m_y_coord = m_falseNorthing - PI * m_radius * sqrt(1.0 - con * con - 2.0 * al * con);
}

void VanDerGrinten::_inverse(double x, double y) 
{

	double xx,yy,xys,c1,c2,c3;
	double a1;
	double m1;
	double con;
	double th1;
	double d;

	/* inverse equations
	-----------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;
	con = PI * m_radius;
	xx = x / con;
	yy = y / con;
	xys = xx * xx + yy * yy;
	c1 = -fabs(yy) * (1.0 + xys);
	c2 = c1 - 2.0 * yy * yy + xx * xx;
	c3 = -2.0 * c1 + 1.0 + 2.0 * yy * yy + xys * xys;
	d = yy * yy / c3 + (2.0 * c2 * c2 * c2 / c3 / c3 / c3 - 9.0 * c1 * c2 / c3 /c3)
		/ 27.0;
	a1 = (c1 - c2 * c2 / 3.0 / c3) / c3;
	m1 = 2.0 * sqrt( -a1 / 3.0);
	con = ((3.0 * d) / a1) / m1;
	if (fabs(con) > 1.0)
	{
	if (con >= 0.0)
		con = 1.0;
	else
		con = -1.0;
	}
	th1 = acos(con) / 3.0;
	if (y >= 0)
	m_latitude = (-m1 * cos(th1 + PI / 3.0) - c2 / 3.0 / c3) * PI;
	else
	m_latitude = -(-m1 * cos(th1 + PI / 3.0) - c2 / 3.0 / c3) * PI;

	if (fabs(xx) < EPSLN)
	{
		m_longitude = m_centerLon;
		return;
	}
	
	m_longitude = Util::adjust_lon(m_centerLon + PI * (xys - 1.0 + sqrt(1.0 + 2.0 * 
								   (xx * xx - yy * yy) + xys * xys)) / 2.0 / xx);
}





