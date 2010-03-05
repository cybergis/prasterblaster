
#include "goodeh.h"

GoodeH::GoodeH() : Projection()
{
	setNumber(GOOD);
	setName("Goode's Homolosine");

	m_centerLons[0] = -1.74532925199;		/* -100.0 degrees */
	m_centerLons[1] = -1.74532925199;		/* -100.0 degrees */
	m_centerLons[2] =  0.523598775598;	/*   30.0 degrees */
	m_centerLons[3] =  0.523598775598;	/*   30.0 degrees */
	m_centerLons[4] = -2.79252680319;		/* -160.0 degrees */
	m_centerLons[5] = -1.0471975512;		/*  -60.0 degrees */
	m_centerLons[6] = -2.79252680319;		/* -160.0 degrees */
	m_centerLons[7] = -1.0471975512;		/*  -60.0 degrees */
	m_centerLons[8] =  0.349065850399;	/*   20.0 degrees */
	m_centerLons[9] =  2.44346095279;		/*  140.0 degrees */
	m_centerLons[10] = 0.349065850399;	/*   20.0 degrees */
	m_centerLons[11] = 2.44346095279;		/*  140.0 degrees */
}

GoodeH::GoodeH(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat)
{
	setNumber(GOOD);
	setName("Goode's Homolosine");

	m_centerLons[0] = -1.74532925199;		/* -100.0 degrees */
	m_centerLons[1] = -1.74532925199;		/* -100.0 degrees */
	m_centerLons[2] =  0.523598775598;	/*   30.0 degrees */
	m_centerLons[3] =  0.523598775598;	/*   30.0 degrees */
	m_centerLons[4] = -2.79252680319;		/* -160.0 degrees */
	m_centerLons[5] = -1.0471975512;		/*  -60.0 degrees */
	m_centerLons[6] = -2.79252680319;		/* -160.0 degrees */
	m_centerLons[7] = -1.0471975512;		/*  -60.0 degrees */
	m_centerLons[8] =  0.349065850399;	/*   20.0 degrees */
	m_centerLons[9] =  2.44346095279;		/*  140.0 degrees */
	m_centerLons[10] = 0.349065850399;	/*   20.0 degrees */
	m_centerLons[11] = 2.44346095279;		/*  140.0 degrees */
}

void GoodeH::_init() 
{
	m_falseEastings[0] = m_radius * -1.74532925199;
	m_falseEastings[1] = m_radius * -1.74532925199;
	m_falseEastings[2] = m_radius * 0.523598775598;
	m_falseEastings[3] = m_radius * 0.523598775598;
	m_falseEastings[4] = m_radius * -2.79252680319;
	m_falseEastings[5] = m_radius * -1.0471975512;
	m_falseEastings[6] = m_radius * -2.79252680319;
	m_falseEastings[7] = m_radius * -1.0471975512;
	m_falseEastings[8] = m_radius * 0.349065850399;
	m_falseEastings[9] = m_radius * 2.44346095279;
	m_falseEastings[10] = m_radius * 0.349065850399;
	m_falseEastings[11] = m_radius * 2.44346095279;
}

void GoodeH::_forward(double lon, double lat)
{
	double delta_lon;	/* Delta longitude (Given longitude - center */
	double theta;
	double delta_theta;
	double constant;
	long i;
	long region;

	/* Forward equations
	-----------------*/
	if (lat >= 0.710987989993)	             /* if on or above 40 44' 11.8" */
	{
		if (lon <= -0.698131700798) 
			region = 0;   /* If to the left of -40 */
		else 
			region = 2;
	}
	else if (lat >= 0.0)			     /* Between 0.0 and 40 44' 11.8" */
	{
		if (lon <= -0.698131700798) 
			region = 1;   /* If to the left of -40 */
		else 
			region = 3;
	}
	else if (lat >= -0.710987989993)   	     /* Between 0.0 & -40 44' 11.8" */
	{
		if (lon <= -1.74532925199) 
			region = 4;  	/* If between -180 and -100 */
		else if (lon <= -0.349065850399) 
			region = 5;	/* If between -100 and -20 */
		else if (lon <= 1.3962634016) 
			region = 8;	/* If between -20 and 80 */
		else 
			region = 9;				/* If between 80 and 180 */
	}
	else						/* Below -40 44' */
	{
		if (lon <= -1.74532925199)
			region = 6;       /* If between -180 and -100 */
		else if (lon <= -0.349065850399) 
			region = 7;     /* If between -100 and -20 */
		else if (lon <= 1.3962634016) 
			region = 10;   /* If between -20 and 80 */
		else 
			region = 11;                            /* If between 80 and 180 */
	}

	if (region==1||region==3||region==4||region==5||region==8||region==9)
	{
		delta_lon = Util::adjust_lon(lon - m_centerLons[region]);
		m_x_coord = m_falseEastings[region] + m_radius * delta_lon * cos(lat);
		m_y_coord = m_radius * lat;
	}
	else
	{
		delta_lon = Util::adjust_lon(lon - m_centerLons[region]);
		theta = lat;
		constant = PI * sin(lat);

		/* Iterate using the Newton-Raphson method to find theta
		-----------------------------------------------------*/
		for (i=0;;i++)
		{
			delta_theta = -(theta + sin(theta) - constant) / (1.0 + cos(theta));
			theta += delta_theta;
			if (fabs(delta_theta) < EPSLN) 
				break;

			if (i >= 50) 
			{
				setError(251);
				return;
			}
		}
		theta /= 2.0;

		/* If the latitude is 90 deg, force the x coordinate to be
			"0 + false easting" this is done here because of precision problems
			with "cos(theta)"
			------------------------------------------------------------------*/
		if (PI / 2 - fabs(lat) < EPSLN)
			delta_lon = 0;
		m_x_coord = m_falseEastings[region] + 0.900316316158 * m_radius * delta_lon * cos(theta);
		m_y_coord = m_radius * (1.4142135623731 * sin(theta) - 0.0528035274542 * Util::sign(lat));
	}
}

void GoodeH::_inverse(double x, double y)
{
	double arg;
	double theta;
	double temp;
	long region;

	/* Inverse equations
	-----------------*/
	if (y >= m_radius * 0.710987989993)                 /* if on or above 40 44' 11.8" */
	{
		if (x <= m_radius * -0.698131700798) 
			region = 0; /* If to the left of -40 */
		else 
			region = 2;
	}
	else if (y >= 0.0)                           /* Between 0.0 and 40 44' 11.8" */
	{
		if (x <= m_radius * -0.698131700798) 
			region = 1; /* If to the left of -40 */
		else 
			region = 3;
	}
	else if (y >= m_radius * -0.710987989993)           /* Between 0.0 & -40 44' 11.8" */
	{
		if (x <= m_radius * -1.74532925199) 
			region = 4;     /* If between -180 and -100 */
		else if (x <= m_radius * -0.349065850399) 
			region = 5; /* If between -100 and -20 */
		else if (x <= m_radius * 1.3962634016) 
			region = 8;  /* If between -20 and 80 */
		else 
			region = 9;                             /* If between 80 and 180 */
	}
	else                                            /* Below -40 44' 11.8" */
	{
		if (x <= m_radius * -1.74532925199) 
			region = 6;     /* If between -180 and -100 */
		else if (x <= m_radius * -0.349065850399) 
			region = 7; /* If between -100 and -20 */
		else if (x <= m_radius * 1.3962634016) 
			region = 10; /* If between -20 and 80 */
		else 
			region = 11;                            /* If between 80 and 180 */
	}
	x = x -  m_falseEastings[region];

	if (region==1||region==3||region==4||region==5||region==8||region==9)
	{
		m_latitude = y / m_radius;
		if (fabs(m_latitude) > HALF_PI) 
		{
			setError(252);
			return;
		}
		temp = fabs(m_latitude) - HALF_PI;
		if (fabs(temp) > EPSLN)
		{
			temp = m_centerLons[region] + x / (m_radius * cos(m_latitude));
			m_longitude = Util::adjust_lon(temp);
		}
		else 
			m_longitude = m_centerLons[region];
	}
	else
	{
		arg = (y + 0.0528035274542 * m_radius * Util::sign(y)) /  (1.4142135623731 * m_radius);
		if (fabs(arg) > 1.0) {
			setError(IN_BREAK);
			return;
		}
		theta = asin(arg);
		m_longitude = Util::adjust_lon(m_centerLons[region]+(x/(0.900316316158 * m_radius * cos(theta))));
		if(m_longitude < -(PI + EPSLN)) {
			setError(IN_BREAK);
			return;
		}
		arg = (2.0 * theta + sin(2.0 * theta)) / PI;
		if (fabs(arg) > 1.0) {
			setError(IN_BREAK);
			return;
		}
		m_latitude = asin(arg);
	}
	/* because of precision problems, long values of 180 deg and -180 deg
	may be mixed.
	----------------------------------------------------------------*/
	if (((x < 0) && (PI - m_longitude < EPSLN)) || ((x > 0) && (PI + m_longitude < EPSLN)))
		m_longitude = -(m_longitude);

	/* Are we in a interrupted area?  If so, return status code of IN_BREAK.
	---------------------------------------------------------------------*/
	if (region == 0 && (m_longitude < -(PI + EPSLN) || m_longitude > -0.698131700798))
		setError(IN_BREAK);
	if (region == 1 && (m_longitude < -(PI + EPSLN) || m_longitude > -0.698131700798))
		setError(IN_BREAK);
	if (region == 2 && (m_longitude < -0.698131700798 || m_longitude > PI + EPSLN))
		setError(IN_BREAK);
	if (region == 3 && (m_longitude < -0.698131700798 || m_longitude > PI + EPSLN))
		setError(IN_BREAK);
	if (region == 4 && (m_longitude < -(PI + EPSLN) || m_longitude > -1.74532925199))
		setError(IN_BREAK);
	if (region == 5 && (m_longitude < -1.74532925199 || m_longitude > -0.349065850399))
		setError(IN_BREAK);
	if (region == 6 && (m_longitude < -(PI + EPSLN) || m_longitude > -1.74532925199))
		setError(IN_BREAK);
	if (region == 7 && (m_longitude < -1.74532925199 || m_longitude > -0.349065850399))
		setError(IN_BREAK);
	if (region == 8 && (m_longitude < -0.349065850399 || m_longitude > 1.3962634016))
		setError(IN_BREAK);
	if (region == 9 && (m_longitude < 1.3962634016|| m_longitude > PI + EPSLN))
		setError(IN_BREAK);
	if (region ==10 && (m_longitude < -0.349065850399 || m_longitude > 1.3962634016))
		setError(IN_BREAK);
	if (region ==11 && (m_longitude < 1.3962634016 || m_longitude > PI + EPSLN))
		setError(IN_BREAK);
}
