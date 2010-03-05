
#include "intmollweide.h"

IntMollweide::IntMollweide() : Projection()
{
	setNumber(IMOLL);
	setName("Interrupted Mollweide");
	m_centerLons[0] = 1.0471975512;
	m_centerLons[1] = -2.96705972839;
	m_centerLons[2] = -0.523598776; 
	m_centerLons[3] = 1.57079632679;
	m_centerLons[4] = -2.44346095279; 
	m_centerLons[5] = -0.34906585;

}

IntMollweide::IntMollweide(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat)
{
	setNumber(IMOLL);
	setName("Interrupted Mollweide");
	m_centerLons[0] = 1.0471975512;
	m_centerLons[1] = -2.96705972839;
	m_centerLons[2] = -0.523598776; 
	m_centerLons[3] = 1.57079632679;
	m_centerLons[4] = -2.44346095279; 
	m_centerLons[5] = -0.34906585;

}

void IntMollweide::_init()
{	
	/* Initialize false eastings for each of the 6 regions
	   ---------------------------------------------------*/
	m_falseEastings[0] = m_radius * -2.19988776387;
	m_falseEastings[1] = m_radius * -0.15713484;
	m_falseEastings[2] = m_radius * 2.04275292359;
	m_falseEastings[3] = m_radius * -1.72848324304;
	m_falseEastings[4] = m_radius * 0.31426968;
	m_falseEastings[5] = m_radius * 2.19988776387;

}



void IntMollweide::_forward(double lon, double lat)
{
	double delta_lon;	/* Delta longitude (Given longitude - center */
	double theta;
	double delta_theta;
	double con;
	long i;
	long region;

	/* Forward equations
	-----------------*/
	/* Note:  PI has been adjusted so that the correct region will be assigned
	when lon = 180 deg.
	----------------------------------------------------------------------*/
	if (lat >= 0.0)
	{
		if (lon >= 0.34906585 && lon < 1.91986217719) 
			region = 0; 
		else if 
			((lon >= 1.919862177 && lon <= (PI + 1.0E-14)) ||
						(lon >= (-PI - 1.0E-14) && lon < -1.745329252))
			region=1; 
		else 
			region = 2;
	}
	else
	{
		if (lon >= 0.34906585 && lon < 2.44346095279) 
			region = 3; 
		else if 
			((lon >= 2.44346095279 && lon <= (PI +1.0E-14)) ||
						(lon >= (-PI - 1.0E-14) && lon<-1.2217304764))
				region=4; 
		else 
			region = 5;
	}

	delta_lon = Util::adjust_lon(lon - m_centerLons[region]);
	theta = lat;
	con = PI * sin(lat);

	/* Iterate using the Newton-Raphson method to find theta
	-----------------------------------------------------*/
	for (i=0;;i++)
	{
		delta_theta = -(theta + sin(theta) - con) / (1.0 + cos(theta));
		theta += delta_theta;
		if (fabs(delta_theta) < EPSLN) 
			break;
		if (i >= 50) {
			setError(2);
			return;
		}
	}
	
	theta /= 2.0;

	/* If the latitude is 90 deg, force the x coordinate to be "0 + false easting"
	this is done here because of percision problems with "cos(theta)"
	--------------------------------------------------------------------------*/
	if (PI / 2 - fabs(lat) < EPSLN)
		delta_lon = 0;

	m_x_coord = m_falseEastings[region] + 0.900316316158 * m_radius * delta_lon * cos(theta);
	m_y_coord = m_radius * 1.4142135623731 * sin(theta);

}

void IntMollweide::_inverse(double x, double y)
{
	
	double theta;
	long region;

	/* Inverse equations
	-----------------*/
	if (y >= 0.0)
	{
		if (x <= m_radius * -1.41421356248) 
			region = 0;
		else if (x <= m_radius * 0.942809042) 
			region = 1;
		else 
			region = 2;
	}
	else
	{
		if (x <= m_radius * -0.942809042) 
			region = 3;
		else if (x <= m_radius * 1.41421356248) 
			region = 4;
		else 
			region = 5;
	}

	x = x - m_falseEastings[region];

	theta = asin(y / (1.4142135623731 * m_radius));
	m_longitude = Util::adjust_lon(m_centerLons[region] + (x / (0.900316316158*m_radius * cos(theta))));
	m_latitude = asin((2.0 * theta + sin(2.0 * theta)) / PI);

	/* Are we in a interrupted area?  If so, return status code of IN_BREAK.
	---------------------------------------------------------------------*/
	if (region == 0 && (m_longitude < 0.34906585 || m_longitude > 1.91986217719))
		setError(IN_BREAK);
	
	if (region == 1 && ((m_longitude < 1.91986217719 && m_longitude > 0.34906585) || 
				(m_longitude > -1.74532925199 && m_longitude < 0.34906585))) 
				setError(IN_BREAK);
	
	if (region == 2 && (m_longitude < -1.745329252 || m_longitude > 0.34906585)) 
		setError(IN_BREAK);
	
	if (region == 3 && (m_longitude < 0.34906585 || m_longitude > 2.44346095279))
		setError(IN_BREAK);
	
	if (region == 4 && ((m_longitude < 2.44346095279 && m_longitude > 0.34906585) || 
				(m_longitude > -1.2217304764 && m_longitude < 0.34906585))) 
				setError(IN_BREAK);

	if (region == 5 && (m_longitude < -1.2217304764 || m_longitude> 0.34906585))
		setError(IN_BREAK);

}
