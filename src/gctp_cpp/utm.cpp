
#include "utm.h"

UTM::UTM() : Projection(), m_scaleFactor(0.0), 
m_e0(0.0), m_e1(0.0), m_e2(0.0), m_e3(0.0), m_e(0.0), m_es(0.0), m_esp(0.0),
m_ml0(0.0), m_ind(0), m_zone(1)
{
	setNumber(_UTM);
	setName("UTM");
	m_scaleFactor = 0.9996;
	m_falseEasting = 500000.0;
	m_centerLat = 0.0;

}

UTM::UTM(double gctpParams[], ProjUnit units, ProjDatum dat, int zone) : 
Projection(gctpParams, units, dat),
m_scaleFactor(0.0), 
m_e0(0.0), m_e1(0.0), m_e2(0.0), m_e3(0.0), m_e(0.0), m_es(0.0), m_esp(0.0),
m_ml0(0.0), m_ind(0), m_zone(zone)
{
	setNumber(_UTM);
	setName("UTM");
	setParamLoad();
	m_scaleFactor = 0.9996;
	m_falseEasting = 500000.0;
	m_centerLat = 0.0;

}


void UTM::_loadFromParams()
{
	Projection::_loadFromParams();
	long err = 0;
	double lon = m_gctpParams[0];
	double lat = m_gctpParams[1];
	
	if(m_zone == 0)
	{

		err = Util::DMSToRad(lon);
		if(err != 0) 
		{
			setError(err);
			return;
		}

		err = Util::DMSToRad(lat);
		if(err != 0)
		{
			setError(err);
			return;
		}

		m_zone = Util::calc_utm_zone(lon * R2D);
		if(lat < 0)
			m_zone *= -1;
	}
}

void UTM::_init() 
{
	double temp;			/* temprorary variables			*/

	if ((abs(m_zone) < 1) || (abs(m_zone) > 60))
	{
		setError(11);
		return;
	}

	m_centerLon = ((6 * abs(m_zone)) - 183) * D2R;
	m_falseNorthing = (m_zone < 0) ? 10000000.0 : 0.0;

	temp = m_rMinor / m_rMajor;
	m_es = 1.0 - SQUARE(temp);
	m_e = sqrt(m_es);
	m_e0 = Util::e0fn(m_es);
	m_e1 = Util::e1fn(m_es);
	m_e2 = Util::e2fn(m_es);
	m_e3 = Util::e3fn(m_es);
	m_ml0 = m_rMajor * Util::mlfn(m_e0, m_e1, m_e2, m_e3, m_centerLat);
	m_esp = m_es / (1.0 - m_es);

	if (m_es < .00001)
		m_ind = 1;
	else 
	    m_ind = 0;

}

void UTM::_forward(double lon, double lat)
{
	double delta_lon;	/* Delta longitude (Given longitude - center 	*/
	double sin_phi, cos_phi;/* sin and cos value				*/
	double al, als;		/* temporary values				*/
	double b;		/* temporary values				*/
	double c, t, tq;	/* temporary values				*/
	double con, n, ml;	/* cone constant, small m			*/

	/* Forward equations
	-----------------*/
	delta_lon = Util::adjust_lon(lon - m_centerLon);
	Util::gctp_sincos(lat, &sin_phi, &cos_phi);

	/* This part was in the fortran code and is for the spherical form 
	----------------------------------------------------------------*/
	if (m_ind != 0)
	{
		b = cos_phi * sin(delta_lon);
		if ((fabs(fabs(b) - 1.0)) < .0000000001)
		{
			setError(93);
			return;
		}
		else
		{
			m_x_coord = .5 * m_rMajor * m_scaleFactor * log((1.0 + b)/(1.0 - b));
			con = acos(cos_phi * cos(delta_lon)/sqrt(1.0 - b*b));
			if (lat < 0)
				con = - con;
			m_y_coord = m_rMajor * m_scaleFactor * (con - m_centerLat); 
			return;
		}
	}

	al  = cos_phi * delta_lon;
	als = SQUARE(al);
	c   = m_esp * SQUARE(cos_phi);
	tq  = tan(lat);
	t   = SQUARE(tq);
	con = 1.0 - m_es * SQUARE(sin_phi);
	n   = m_rMajor / sqrt(con);
	ml  = m_rMajor * Util::mlfn(m_e0, m_e1, m_e2, m_e3, lat);

	m_x_coord  = m_scaleFactor * n * al * (1.0 + als / 6.0 * (1.0 - t + c + als / 20.0 *
		(5.0 - 18.0 * t + SQUARE(t) + 72.0 * c - 58.0 * m_esp))) + m_falseEasting;

	m_y_coord  = m_scaleFactor * (ml - m_ml0 + n * tq * (als * (0.5 + als / 24.0 *
		(5.0 - t + 9.0 * c + 4.0 * SQUARE(c) + als / 30.0 * (61.0 - 58.0 * t
		+ SQUARE(t) + 600.0 * c - 330.0 * m_esp))))) + m_falseNorthing;
}

void UTM::_inverse(double x, double y)
{
	double con,phi;		/* temporary angles				*/
	double delta_phi;	/* difference between longitudes		*/
	long i;			/* counter variable				*/
	double sin_phi, cos_phi, tan_phi;	/* sin cos and tangent values	*/
	double c, cs, t, ts, n, r, d, ds;	/* temporary variables		*/
	double f, h, g, temp;			/* temporary variables		*/
	long max_iter = 6;			/* maximun number of iterations	*/

	/* fortran code for spherical form 
	--------------------------------*/
	if (m_ind != 0)
	{
		f = exp(x/(m_rMajor * m_scaleFactor));
		g = .5 * (f - 1/f);
		temp = m_centerLat + y/(m_rMajor * m_scaleFactor);
		h = cos(temp);
		con = sqrt((1.0 - h * h)/(1.0 + g * g));
		m_latitude = Util::asinz(con);
		if (temp < 0)
			m_latitude = -m_latitude;
		if ((g == 0) && (h == 0))
		{
			m_longitude = m_centerLon;
			return;
		}
		else
		{
			m_longitude = Util::adjust_lon(atan2(g,h) + m_centerLon);
			return;
		}
	}

	/* Inverse equations
	-----------------*/
	x = x - m_falseEasting;
	y = y - m_falseNorthing;

	con = (m_ml0 + y / m_scaleFactor) / m_rMajor;
	phi = con;
	for (i=0;;i++)
	{
		delta_phi=((con + m_e1 * sin(2.0*phi) - m_e2 * sin(4.0*phi) + m_e3 * sin(6.0*phi))
			/ m_e0) - phi;
	
		phi += delta_phi;
		if (fabs(delta_phi) <= EPSLN) 
			break;

		if (i >= max_iter) 
		{ 
			setError(95);
			return;
		}
	}
	if (fabs(phi) < HALF_PI)
	{
		Util::gctp_sincos(phi, &sin_phi, &cos_phi);
		tan_phi = tan(phi);
		c    = m_esp * SQUARE(cos_phi);
		cs   = SQUARE(c);
		t    = SQUARE(tan_phi);
		ts   = SQUARE(t);
		con  = 1.0 - m_es * SQUARE(sin_phi); 
		n    = m_rMajor / sqrt(con);
		r    = n * (1.0 - m_es) / con;
		d    = x / (n * m_scaleFactor);
		ds   = SQUARE(d);

		m_latitude = phi - (n * tan_phi * ds / r) * (0.5 - ds / 24.0 * (5.0 + 3.0 * t + 
			10.0 * c - 4.0 * cs - 9.0 * m_esp - ds / 30.0 * (61.0 + 90.0 * t +
			298.0 * c + 45.0 * ts - 252.0 * m_esp - 3.0 * cs)));
		m_longitude = Util::adjust_lon(m_centerLon + (d * (1.0 - ds / 6.0 * (1.0 + 2.0 * t +
			c - ds / 20.0 * (5.0 - 2.0 * c + 28.0 * t - 3.0 * cs + 8.0 * m_esp +
			24.0 * ts))) / cos_phi));
	}
	else
	{
		m_latitude = HALF_PI * Util::sign(y);
		m_longitude = m_centerLon;
	}
}

