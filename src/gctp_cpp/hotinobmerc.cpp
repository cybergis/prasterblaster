
#include "hotinobmerc.h"

HotineObMerc::HotineObMerc(): Projection(), m_azimuth(0.0), m_scaleFactor(1),
m_e(0.0), m_es(0.0), m_sinCenterLat(0.0), m_cosCenterLat(0.0), m_bl(0.0),
m_al(0.0), m_ts(0.0), m_d(0.0), m_el(0.0), m_u(0.0), m_singam(0.0), m_cosgam(0.0),
m_sinaz(0.0), m_cosaz(0.0), m_lat1(0.0), m_lat2(0.0), m_lon1(0.0), m_lon2(0.0),
m_mode(0)
{
	setNumber(HOM);
	setName("Hotine Oblique Mercator");
}

HotineObMerc::HotineObMerc(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat), m_azimuth(0.0), m_scaleFactor(1),
m_e(0.0), m_es(0.0), m_sinCenterLat(0.0), m_cosCenterLat(0.0), m_bl(0.0),
m_al(0.0), m_ts(0.0), m_d(0.0), m_el(0.0), m_u(0.0), m_singam(0.0), m_cosgam(0.0),
m_sinaz(0.0), m_cosaz(0.0), m_lat1(0.0), m_lat2(0.0), m_lon1(0.0), m_lon2(0.0),
m_mode(0)
{
	setNumber(HOM);
	setName("Hotine Oblique Mercator");
	setParamLoad();
}

void HotineObMerc::_init()
{
	double temp;			/* temporary variable		*/
	double con,com;
	double h,l,ts1,ts2;
	double j,p,dlon;
	double f,g,gama;
	double sinphi;


	temp = m_rMinor / m_rMajor;
	m_es = 1.0 - SQUARE(temp);
	m_e = sqrt(m_es);

	Util::gctp_sincos(m_centerLat,&m_cosCenterLat,&m_cosCenterLat);
	con = 1.0 - m_es * m_cosCenterLat * m_cosCenterLat;
	com = sqrt(1.0 - m_es);
	m_bl = sqrt(1.0 + m_es * pow(m_cosCenterLat,4.0)/(1.0 - m_es));
	m_al = m_rMajor * m_bl * m_scaleFactor * com / con;
	if (fabs(m_centerLat) < EPSLN)
	{
		m_ts = 1.0;
		m_d = 1.0;
		m_el = 1.0;
	}
	else
	{
		m_ts = Util::tsfnz(m_e,m_centerLat,m_cosCenterLat);
		con = sqrt(con);
		m_d = m_bl * com / (m_cosCenterLat * con);
		if ((m_d * m_d - 1.0) > 0.0)
		{
			if (m_centerLat >= 0.0)
				f = m_d + sqrt(m_d * m_d - 1.0);
			else
				f = m_d - sqrt(m_d * m_d - 1.0);
		}
		else
			f = m_d;
		m_el = f * pow(m_ts,m_bl);
	}


	if (m_mode != 0)
	{
		g = .5 * (f - 1.0/f);
		gama = Util::asinz(sin(m_azimuth) / m_d);
		m_centerLon = m_centerLon - Util::asinz(g * tan(gama))/m_bl;

		con = fabs(m_centerLat);
		if ((con > EPSLN) && (fabs(con - HALF_PI) > EPSLN))
		{
			Util::gctp_sincos(gama,&m_singam,&m_cosgam);
			Util::gctp_sincos(m_azimuth,&m_sinaz,&m_cosaz);
			if (m_centerLat >= 0)
				m_u =  (m_al / m_bl) * atan(sqrt(m_d*m_d - 1.0)/m_cosaz);
			else
				m_u =  -(m_al / m_bl) * atan(sqrt(m_d*m_d - 1.0)/m_cosaz);
		}
		else
		{
			setError(201);
			return;
		}
	}
	else
	{
		sinphi = sin(m_lat1);
		ts1 = Util::tsfnz(m_e,m_lat1,sinphi);
		sinphi = sin(m_lat2);
		ts2 = Util::tsfnz(m_e,m_lat2,sinphi);
		h = pow(ts1,m_bl);
		l = pow(ts2,m_bl);
		f = m_el/h;
		g = .5 * (f - 1.0/f);
		j = (m_el * m_el - l * h)/(m_el * m_el + l * h);
		p = (l - h) / (l + h);
		dlon = m_lon1 - m_lon2;
		if (dlon < -PI)
			m_lon2 = m_lon2 - 2.0 * PI;
		if (dlon > PI)
			m_lon2 = m_lon2 + 2.0 * PI;
		dlon = m_lon1 - m_lon2;
		m_centerLon = .5 * (m_lon1 + m_lon2) - atan(j * tan(.5 * m_bl * dlon)/p)/m_bl;
		dlon  = Util::adjust_lon(m_lon1 - m_centerLon);
		gama = atan(sin(m_bl * dlon)/g);
		m_azimuth = Util::asinz(m_d * sin(gama));

		if (fabs(m_lat1 - m_lat2) <= EPSLN)
		{
			setError(202);
			return;
		}
		else
			con = fabs(m_lat1);
		if ((con <= EPSLN) || (fabs(con - HALF_PI) <= EPSLN))
		{
			setError(202);
			return;
		}
		else 
			if (fabs(fabs(m_centerLat) - HALF_PI) <= EPSLN)
			{
				setError(202);
				return;
			}

			Util::gctp_sincos(gama,&m_singam,&m_cosgam);
			Util::gctp_sincos(m_azimuth,&m_sinaz,&m_cosaz);
			if (m_centerLat >= 0)
				m_u =  (m_al/m_bl) * atan(sqrt(m_d * m_d - 1.0)/m_cosaz);
			else
				m_u = -(m_al/m_bl) * atan(sqrt(m_d * m_d - 1.0)/m_cosaz);
	}

}

void HotineObMerc::_forward(double lon, double lat)
{
	double sin_phi;/* sin and cos value				*/
	double t;	/* temporary values				*/
	double con;	/* cone constant, small m			*/
	double q,us,vl;
	double ul,vs;
	double s;
	double dlon;
	double ts1;

	/* Forward equations
	-----------------*/
	sin_phi = sin(lat);
	dlon = Util::adjust_lon(lon - m_centerLon);
	vl = sin(m_bl * dlon);
	if (fabs(fabs(lat) - HALF_PI) > EPSLN)
	{
		ts1 = Util::tsfnz(m_e,lat,sin_phi);
		q = m_el / (pow(ts1,m_bl));
		s = .5 * (q - 1.0 / q);
		t = .5 * (q + 1.0/ q);
		ul = (s * m_singam - vl * m_cosgam) / t;
		con = cos(m_bl * dlon);
		if (fabs(con) < .0000001)
		{
			us = m_al * m_bl * dlon;
		}
		else
		{
			us = m_al * atan((s * m_cosgam + vl * m_singam) / con)/m_bl;
			if (con < 0)
				us = us + PI * m_al / m_bl;
		}
	}
	else
	{
		if (lat >= 0)
			ul = m_singam;
		else
			ul = -m_singam;

		us = m_al * lat / m_bl;
	}
	if (fabs(fabs(ul) - 1.0) <= EPSLN)
	{
		setError(205);
		return;
	}

	vs = .5 * m_al * log((1.0 - ul)/(1.0 + ul)) / m_bl;
	us = us - m_u;
	m_x_coord = m_falseEasting + vs * m_cosaz + us * m_sinaz;
	m_y_coord = m_falseNorthing + us * m_cosaz - vs * m_sinaz;
}

void HotineObMerc::_inverse(double x, double y)
{
	double theta;		/* angle					*/
	double t;	/* temporary values				*/
	double con;	/* cone constant, small m			*/
	double vs,us,q,s,ts1;
	double vl,ul;
	long   flag;

	/* Inverse equations
	-----------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;
	flag = 0;
	vs = x * m_cosaz - y * m_sinaz;
	us = y * m_cosaz + x * m_sinaz;
	us = us + m_u;
	q = exp(-m_bl * vs / m_al);
	s = .5 * (q - 1.0/q);
	t = .5 * (q + 1.0/q);
	vl = sin(m_bl * us / m_al);
	ul = (vl * m_cosgam + s * m_singam)/t;
	if (fabs(fabs(ul) - 1.0) <= EPSLN)
	{
		m_longitude = m_centerLon;
		if (ul >= 0.0)
			m_latitude = HALF_PI;
		else
			m_latitude = -HALF_PI;
	}
	else
	{
		con = 1.0 / m_bl;
		ts1 = pow((m_el / sqrt((1.0 + ul) / (1.0 - ul))),con);
		m_latitude = Util::phi2z(m_e,ts1,&flag);
		if (flag != 0) 
		{
			setError(flag);
			return;
		}

		con = cos(m_bl * us /m_al);
		theta = m_centerLon - atan2((s * m_cosgam - vl * m_singam) , con)/m_bl;
		m_longitude = Util::adjust_lon(theta);
	}
}

void HotineObMerc::_loadFromParams()
{
	Projection::_loadFromParams();
	setScaleFactor(m_gctpParams[2]);
	setMode((int)m_gctpParams[12]);
	if(m_mode == 1) 
		setAzimuth(m_gctpParams[3]);
	
	else
	{
		setLon1(m_gctpParams[8]);
		setLat1(m_gctpParams[9]);
		setLon2(m_gctpParams[10]);
		setLat2(m_gctpParams[11]);
	}
}

void HotineObMerc::setLon1(double lon) 
{
	convertAndSetAngle(m_lon1, lon);
	setInit();

}

void HotineObMerc::setLat1(double lat) 
{
	convertAndSetAngle(m_lat1, lat);
	setInit();
}

void HotineObMerc::setLat2(double lat) 
{
	convertAndSetAngle(m_lat2, lat);
	setInit();
}

void HotineObMerc::setLon2(double lon) 
{
	convertAndSetAngle(m_lon2, lon);
	setInit();
}

void HotineObMerc::setAzimuth(double angle)
{
	convertAndSetAngle(m_azimuth, angle);
	setInit();
}

