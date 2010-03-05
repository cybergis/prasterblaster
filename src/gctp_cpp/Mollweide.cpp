#include "mollweide.h"

Mollweide::Mollweide(): Projection() 
{
	setName("Mollweide");
	setNumber(MOLL);
}

Mollweide::Mollweide(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat) 
{
	setName("Mollweide");
	setNumber(MOLL);
}

void Mollweide::_init() 
{
	return;
}

void Mollweide::_forward(double lon, double lat) 
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
  con = PI * sin(lat);

 /* Iterate using the Newton-Raphson method to find theta
  -----------------------------------------------------*/
  for (i=0;;i++)
	{
	delta_theta = -(theta + sin(theta) - con)/ (1.0 + cos(theta));
	theta += delta_theta;
	if (fabs(delta_theta) < EPSLN) break;
	if (i >= 50) {
		setError(241);
		return;
	}
	  
	}
	theta /= 2.0;

	/*	If the latitude is 90 deg, force the x coordinate to be "0 + false easting"
	this is done here because of precision problems with "cos(theta)"
	 --------------------------------------------------------------------------*/
	if (PI/2 - fabs(lat) < EPSLN)
		delta_lon =0;
	
	m_x_coord = 0.900316316158 * m_radius * delta_lon * cos(theta) + m_falseEasting;
	m_y_coord = 1.4142135623731 * m_radius * sin(theta) + m_falseNorthing;

}

void Mollweide::_inverse(double x, double y) 
{
  double theta;
  double arg;

  /* Inverse equations
   -----------------*/
  x -= m_falseEasting; 
  y -= m_falseNorthing;
  arg = y /  (1.4142135623731 * m_radius);

  /* Because of division by zero problems, 'arg' can not be 1.0.  Therefore
     a number very close to one is used instead.
     -------------------------------------------------------------------*/
  if(fabs(arg) > 0.999999999999) 
	  arg=0.999999999999;
  theta = asin(arg);
  m_longitude = Util::adjust_lon(m_centerLon + (x / (0.900316316158 * m_radius * cos(theta))));
  if(m_longitude < (-PI)) 
	  m_longitude = -PI;
  if(m_longitude > PI) 
	  m_longitude = PI;
  arg = (2.0 * theta + sin(2.0 * theta)) / PI;
  if(fabs(arg) > 1.0)
	  arg=1.0;
  m_latitude = asin(arg);

}

