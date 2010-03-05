
#include "equidistantc.h"

EquidistantC::EquidistantC(): Projection(), m_e(0.0), m_es(0.0),
m_e0(0.0), m_e1(0.0), m_e2(0.0), m_e3(0.0), m_ml0(0.0), m_ns(0.0),
m_g(0.0), m_rh(0.0)
{
	setNumber(EQUIDC);
	setName("Equidistant Conic");
}

EquidistantC::EquidistantC(double gctpParams[], ProjUnit units, ProjDatum dat): 
Projection(gctpParams, units, dat), m_e(0.0), m_es(0.0), 
m_e0(0.0), m_e1(0.0), m_e2(0.0), m_e3(0.0), m_ml0(0.0), m_ns(0.0),
m_g(0.0), m_rh(0.0)
{
	setNumber(EQUIDC);
	setName("Equidistant Conic");
	setParamLoad();
}

void EquidistantC::_init() 
{
	double temp;			/* temporary variable		*/
	double sinphi,cosphi;		/* sin and cos values		*/
	double ms1,ms2;
	double ml1,ml2;

	temp = m_rMinor / m_rMajor;
	m_es = 1.0 - SQUARE(temp);
	m_e = sqrt(m_es);
	m_e0 = Util::e0fn(m_es);
	m_e1 = Util::e1fn(m_es);
	m_e2 = Util::e2fn(m_es);
	m_e3 = Util::e3fn(m_es);

	Util::gctp_sincos(m_stdParallelLat1,&sinphi,&cosphi);
	ms1 = Util::msfnz(m_e,sinphi,cosphi);
	ml1 = Util::mlfn(m_e0, m_e1, m_e2, m_e3, m_stdParallelLat1);
	
	if (m_mode != 0)
	{
		if (fabs(m_stdParallelLat1 + m_stdParallelLat2) < EPSLN) {
			setError(81);
			return;
		}
		
		Util::gctp_sincos(m_stdParallelLat2,&sinphi,&cosphi);
		ms2 = Util::msfnz(m_e,sinphi,cosphi);
		ml2 = Util::mlfn(m_e0, m_e1, m_e2, m_e3, m_stdParallelLat2);
		if (fabs(m_stdParallelLat1 - m_stdParallelLat2) >= EPSLN)
			m_ns = (ms1 - ms2) / (ml2 - ml1);
		else
			m_ns = sinphi;
	}
	else
		m_ns = sinphi;	

	m_g = ml1 + ms1/m_ns;
	m_ml0 = Util::mlfn(m_e0, m_e1, m_e2, m_e3, m_centerLat);
	m_rh = m_rMajor * (m_g - m_ml0);


}

void EquidistantC::_forward(double lon, double lat)
{
	double ml;
	double theta;
	double rh1;
	
	/* Forward equations
	-----------------*/
	ml = Util::mlfn(m_e0, m_e1, m_e2, m_e3, lat);
	rh1 = m_rMajor * (m_g - ml);
	theta = m_ns * Util::adjust_lon(lon - m_centerLon);
	m_x_coord = m_falseEasting  + rh1 * sin(theta);
	m_y_coord = m_falseNorthing + m_rh - rh1 * cos(theta);

}

void EquidistantC::_inverse(double x, double y)
{
	double rh1;
	double ml;
	double con;
	double theta;
	long   flag;

	flag = 0;
	x -= m_falseEasting;
	y  = m_rh - y + m_falseNorthing;
	if (m_ns >= 0)
	{
		rh1 = sqrt(x * x + y * y); 
		con = 1.0;
	}
	else
	{
		rh1 = -sqrt(x * x + y * y); 
		con = -1.0;
	}

	theta = 0.0;
	if (rh1  != 0.0)
		theta = atan2(con * x, con * y);
	
	ml = m_g - rh1 / m_rMajor;
	m_latitude = Util::phi3z(ml,m_e0,m_e1,m_e2,m_e3,&flag);
	
	if (flag != 0) {
		setError(flag);
		return;
	}

	m_longitude = Util::adjust_lon(m_centerLon + theta / m_ns);

}

void EquidistantC::_loadFromParams() {
	Projection::_loadFromParams();
	setMode((int)m_gctpParams[8]);
}
	

