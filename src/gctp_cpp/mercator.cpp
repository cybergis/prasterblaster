
#include "mercator.h"

Mercator::Mercator(): Projection(), m_e(0.0),
m_es(0.0), m_m1(0.0)
{
	setName("Mercator");
	setNumber(MERCAT);
}

Mercator::Mercator(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat), m_e(0.0),
m_es(0.0), m_m1(0.0)
{
		setName("Mercator");
		setNumber(MERCAT);
}

void Mercator::_init() {
	double temp;			/* temporary variable		*/

	temp = m_rMinor / m_rMajor;
	m_es = 1.0 - SQUARE(temp);
	m_e = sqrt(m_es);
	m_m1 = cos(m_centerLat)/(sqrt(1.0 - m_es * sin(m_centerLat) * sin(m_centerLat)));

}
void Mercator::_inverse(double x, double y)
{
	double ts;		/* small t value				*/
	long flag;		/* error flag 					*/

	/* Inverse equations
	-----------------*/
	flag = 0;
	x -= m_falseEasting;
	y -= m_falseNorthing;
	ts = exp(-y/(m_rMajor * m_m1));
	m_latitude = Util::phi2z(m_e,ts,&flag);
	
	if (flag != 0) {
		setError(flag);
		return;
	}

	m_longitude = Util::adjust_lon(m_centerLon + x/(m_rMajor * m_m1));

}

void Mercator::_forward(double lon, double lat)
{
	double ts;		/* small t value				*/
	double sinphi;		/* sin value					*/


	/* Forward equations
	 -----------------*/
	if (fabs(fabs(lat) - HALF_PI)  <= EPSLN) {
		setError(52);
		return;
	}

	else
	 {
		sinphi = sin(lat);
		ts = Util::tsfnz(m_e,lat,sinphi);
		m_x_coord = m_falseEasting + m_rMajor * m_m1 * Util::adjust_lon(lon - m_centerLon);
		m_y_coord = m_falseNorthing - m_rMajor * m_m1 * log(ts);
	 }

}













