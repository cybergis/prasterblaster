
#include "oblatedeqarea.h"

OblatedEqArea::OblatedEqArea() : Projection(),
m_theta(0.0), m_m(0.0), m_n(0.0), m_sinLatO(0.0), m_cosLatO(0.0)
{
	setNumber(OBEQA);
	setName("Oblated Equal Area");
}

OblatedEqArea::OblatedEqArea(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat),
m_theta(0.0), m_m(0.0), m_n(0.0), m_sinLatO(0.0), m_cosLatO(0.0)
{
	setNumber(OBEQA);
	setName("Oblated Equal Area");
	setParamLoad();
}

void OblatedEqArea::_loadFromParams() 
{
	Projection::_loadFromParams();
	setShapeM(m_gctpParams[2]);
	setShapeN(m_gctpParams[3]);
	setAngle(m_gctpParams[8]);
}

void OblatedEqArea::_init()
{
	
	// Calculate the sine and cosine of the latitude of the center of the map
	Util::gctp_sincos(m_centerLat, &m_sinLatO, &m_cosLatO);

}



void OblatedEqArea::setAngle(double theta)
{
	long iflg = 0;
	double temp = Util::paksz(theta,&iflg)* 3600 * S2R;
	if(iflg != 0) 
	{
		setError(iflg);
		return;
	}
	m_theta = temp;
}

void OblatedEqArea::_forward(double lon, double lat) 
{
	double delta_lon;
	double sin_delta_lon;
	double cos_delta_lon;
	double sin_lat;
	double cos_lat;
	double z;
	double Az;
	double sin_Az;
	double cos_Az;
	double temp;			/* Re-used temporary variable */
	double x_prime;
	double y_prime;
	double M;
	double N;


	/* Forward equations
	-----------------*/
	delta_lon = lon - m_centerLon;
	Util::gctp_sincos(lat, &sin_lat, &cos_lat);
	Util::gctp_sincos(delta_lon, &sin_delta_lon, &cos_delta_lon);
	z = acos(m_sinLatO * sin_lat + m_cosLatO * cos_lat * cos_delta_lon);
	Az = atan2(cos_lat * sin_delta_lon , m_cosLatO * sin_lat - m_sinLatO * 
		cos_lat * cos_delta_lon) + m_theta;
	Util::gctp_sincos(Az, &sin_Az, &cos_Az);
	temp = 2.0 * sin(z / 2.0);
	x_prime = temp * sin_Az;
	y_prime = temp * cos_Az;
	M = asin(x_prime / 2.0);
	temp = y_prime / 2.0 * cos(M) / cos(2.0 * M / m_m);
	N = asin(temp);
	m_y_coord = m_n * m_radius * sin(2.0 * N / m_n) + m_falseEasting;
	m_x_coord = m_m * m_radius * sin(2.0 * M / m_m) * cos(N) / cos(2.0 * N / m_n) + m_falseNorthing;

}

void OblatedEqArea::_inverse(double x, double y)
{
	double z;
	double sin_z;
	double cos_z;
	double Az;
	double temp;			/* Re-used temporary variable */
	double x_prime;
	double y_prime;
	double M;
	double N;
	double diff_angle;
	double sin_diff_angle;
	double cos_diff_angle;


	/* Inverse equations
	-----------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;
	N = (m_n / 2.0) * asin(y / (m_n * m_radius));
	temp = x / (m_m * m_radius) * cos(2.0 * N / m_n) / cos(N);
	M = (m_m / 2.0) * asin(temp);
	x_prime = 2.0 * sin(M);
	y_prime = 2.0 * sin(N) * cos(2.0 * M / m_m) / cos(M);
	temp = sqrt(x_prime * x_prime + y_prime * y_prime) / 2.0;
	z = 2.0 * asin(temp);
	Az = atan2(x_prime, y_prime);
	diff_angle = Az - m_theta;
	Util::gctp_sincos(diff_angle, &sin_diff_angle, &cos_diff_angle);
	Util::gctp_sincos(z, &sin_z, &cos_z);
	m_latitude = asin(m_sinLatO * cos_z + m_cosLatO * sin_z * cos_diff_angle);
	m_longitude = Util::adjust_lon(m_centerLon + atan2((sin_z * sin_diff_angle), (m_cosLatO *
								   cos_z - m_sinLatO * sin_z * cos_diff_angle)));

}