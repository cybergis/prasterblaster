
#include "lambertcc.h"

LambertCC::LambertCC() : Projection(), m_es(0.0), m_e(0.0), m_ns(0.0),
m_f0(0.0), m_rh(0.0)
{
	setNumber(LAMCC);
	setName("Lambert Conformal Conic");
}

LambertCC::LambertCC(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat), m_es(0.0), m_e(0.0), m_ns(0.0),
m_f0(0.0), m_rh(0.0)
{
	setNumber(LAMCC);
	setName("Lambert Conformal Conic");
}

void LambertCC::_init() 
{
	double sin_po;                  /* sin value                            */
	double cos_po;                  /* cos value                            */
	double con;                     /* temporary variable                   */
	double ms1;                     /* small m 1                            */
	double ms2;                     /* small m 2                            */
	double temp;                    /* temporary variable                   */
	double ts0;                     /* small t 0                            */
	double ts1;                     /* small t 1                            */
	double ts2;                     /* small t 2                            */

	/* Standard Parallels cannot be equal and on opposite sides of the equator
	------------------------------------------------------------------------*/
	if (fabs(m_stdParallelLat1+m_stdParallelLat2) < EPSLN)
	{
		setError(41);
		return;
	}
	   
	temp = m_rMinor / m_rMajor;
	m_es = 1.0 - SQUARE(temp);
	m_e = sqrt(m_es);

	Util::gctp_sincos(m_stdParallelLat1,&sin_po,&cos_po);
	con = sin_po;
	ms1 = Util::msfnz(m_e,sin_po,cos_po);
	ts1 = Util::tsfnz(m_e,m_stdParallelLat1,sin_po);
	Util::gctp_sincos(m_stdParallelLat2,&sin_po,&cos_po);
	ms2 = Util::msfnz(m_e,sin_po,cos_po);
	ts2 = Util::tsfnz(m_e,m_stdParallelLat2,sin_po);
	sin_po = sin(m_centerLat);
	ts0 = Util::tsfnz(m_e,m_centerLat,sin_po);

	if (fabs(m_stdParallelLat1 - m_stdParallelLat2) > EPSLN)
		m_ns = log (ms1/ms2)/ log (ts1/ts2);
	else
		m_ns = con;
	m_f0 = ms1 / (m_ns * pow(ts1,m_ns));
	m_rh = m_rMajor * m_f0 * pow(ts0,m_ns);

}

void LambertCC::_forward(double lon, double lat)
{
	double con;                     /* temporary angle variable             */
	double rh1;                     /* height above ellipsoid               */
	double sinphi;                  /* sin value                            */
	double theta;                   /* angle                                */
	double ts;                      /* small value t                        */
	
	con  = fabs( fabs(lat) - HALF_PI);
	if (con > EPSLN)
	{
		sinphi = sin(lat);
		ts = Util::tsfnz(m_e,lat,sinphi);
		rh1 = m_rMajor * m_f0 * pow(ts,m_ns);
	}
	else
	{
		con = lat * m_ns;
		if (con <= 0)
		{
			setError(44);
			return;
		}
		rh1 = 0;
	}

	theta = m_ns * Util::adjust_lon(lon - m_centerLon);
	m_x_coord = rh1 * sin(theta) + m_falseEasting;
	m_y_coord = m_rh - rh1 * cos(theta) + m_falseNorthing;	

}

void LambertCC::_inverse(double x, double y)
{
	double rh1;			/* height above ellipsoid	*/
	double con;			/* sign variable		*/
	double ts;			/* small t			*/
	double theta;			/* angle			*/
	long   flag;			/* error flag			*/

	flag = 0;
	x -= m_falseEasting;
	y = m_rh - y + m_falseNorthing;
	if (m_ns > 0)
		{
		rh1 = sqrt (x * x + y * y);
		con = 1.0;
		}
	else
		{
		rh1 = -sqrt (x * x + y * y);
		con = -1.0;
		}
	theta = 0.0;
	if (rh1 != 0)
		theta = atan2((con * x),(con * y));
	if ((rh1 != 0) || (m_ns > 0.0))
	{
		con = 1.0/m_ns;
		ts = pow((rh1/(m_rMajor * m_f0)),con);
		m_latitude = Util::phi2z(m_e,ts,&flag);
		if (flag != 0)
		{
			setError(flag);
			return;
		}
	}
	else
		m_latitude = -HALF_PI;
	
	m_longitude = Util::adjust_lon(theta/m_ns + m_centerLon);

}
 



