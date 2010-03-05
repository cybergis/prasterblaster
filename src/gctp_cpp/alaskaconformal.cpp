
#include "alaskaconformal.h"

AlaskaConformal::AlaskaConformal() : 
Projection(), m_sinCenterLat(0.0),
m_cosCenterLat(0.0), m_e(0.0), m_n(6)
{
	double es;
	double chi;
	double esphi;

	setNumber(ALASKA);
	setName("Alaska Conformal");
	m_centerLon = -152.0 * D2R;
	m_centerLat = 64.0 * D2R;
	m_acoef[0] = 0;
	m_acoef[1]= 0.9945303;  
	m_acoef[2]= 0.0052083;   
	m_acoef[3]= 0.0072721;    
	m_acoef[4]= -0.0151089;    
	m_acoef[5]= 0.0642675;      
	m_acoef[6]= 0.3582802;
	m_bcoef[0] = 0;
	m_bcoef[1]= 0.0;      
	m_bcoef[2]= -.0027404; 
	m_bcoef[3]= 0.0048181;  
	m_bcoef[4]= -0.1932526;  
	m_bcoef[5]= -0.1381226;
	m_bcoef[6]= -0.2884586; 

	es = .006768657997291094;
	m_e = sqrt(es);
	esphi = m_e * sin(m_centerLat);
	
	chi = 2.0 * atan(tan((HALF_PI + m_centerLat)/2.0) * 
				pow(((1.0 - esphi)/(1.0 + esphi)),(m_e/2.0))) - HALF_PI;
	
	Util::gctp_sincos(chi,&m_sinCenterLat,&m_cosCenterLat);

}

AlaskaConformal::AlaskaConformal(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat), m_sinCenterLat(0.0),
m_cosCenterLat(0.0), m_e(0.0), m_n(6)

{
	double es;
	double chi;
	double esphi;

	setNumber(ALASKA);
	setName("Alaska Conformal");
	m_centerLon = -152.0 * D2R;
	m_centerLat = 64.0 * D2R;

	m_acoef[0] = 0;
	m_acoef[1]= 0.9945303;  
	m_acoef[2]= 0.0052083;   
	m_acoef[3]= 0.0072721;    
	m_acoef[4]= -0.0151089;    
	m_acoef[5]= 0.0642675;      
	m_acoef[6]= 0.3582802;
	m_bcoef[0] = 0;
	m_bcoef[1]= 0.0;      
	m_bcoef[2]= -.0027404; 
	m_bcoef[3]= 0.0048181;  
	m_bcoef[4]= -0.1932526;  
	m_bcoef[5]= -0.1381226;
	m_bcoef[6]= -0.2884586; 

	es = .006768657997291094;
	m_e = sqrt(es);
	esphi = m_e * sin(m_centerLat);
	
	chi = 2.0 * atan(tan((HALF_PI + m_centerLat)/2.0) * 
				pow(((1.0 - esphi)/(1.0 + esphi)),(m_e/2.0))) - HALF_PI;
	
	Util::gctp_sincos(chi,&m_sinCenterLat,&m_cosCenterLat);

}

void AlaskaConformal::_init()
{
	return;
}

void AlaskaConformal::_forward(double lon, double lat)
{

	double dlon;
	double sinlon,coslon;
	double sinphi,cosphi;
	double esphi;
	double g;
	double s;
	double xp;
	double yp;
	double ar;
	double ai;
	double br;
	double bi;
	double arn;
	double ain;
	double chi;
	double r;
	long j;


	/* Forward equations
	-----------------*/
	dlon = Util::adjust_lon( lon - m_centerLon);

	/* caluclate x' and y' for Oblique Stereographic Proj for LAT/LONG
	----------------------------------------------------------------*/
	Util::gctp_sincos(dlon,&sinlon,&coslon);
	esphi = m_e * sin(lat);
	chi = 2.0 * atan(tan((HALF_PI + lat) / 2.0) * 
		pow(((1.0 - esphi) / (1.0 + esphi)),(m_e/2.0))) - HALF_PI;

	Util::gctp_sincos(chi,&sinphi,&cosphi);
	g = m_sinCenterLat * sinphi + m_cosCenterLat * cosphi * coslon;
	s = 2.0 / (1.0 + g);
	xp = s * cosphi * sinlon;
	yp = s * (m_cosCenterLat * sinphi - m_sinCenterLat * cosphi * coslon);

	/* Use Knuth algorithm for summing complex terms, to convert
	Oblique Stereographic to Modified-Stereographic coord
	----------------------------------------------------------*/
	r = xp + xp;
	s = xp*xp + yp*yp;
	ar = m_acoef[m_n];
	ai = m_bcoef[m_n];
	br = m_acoef[m_n -1];
	bi = m_bcoef[m_n -1];
	for (j =2; j <= m_n; j++)
	{
		arn = br + r * ar;
		ain = bi + r * ai; 
		if (j < m_n)
		{
			br = m_acoef[m_n - j] - s * ar;
			bi = m_bcoef[m_n - j] - s * ai;
			ar = arn;
			ai = ain;
		}
	}
	br = -s * ar;
	bi = -s * ai;
	ar = arn;
	ai = ain;
	m_x_coord = (xp * ar - yp * ai + br) * m_rMajor + m_falseEasting;
	m_y_coord = (yp * ar + xp * ai + bi) * m_rMajor + m_falseNorthing;
}

void AlaskaConformal::_inverse(double x, double y)
{
	double esphi;
	double r;
	double s;
	double br;
	double bi;
	double ai;
	double ar;
	double ci;
	double cr;
	double di;
	double dr;
	double arn;
	double ain;
	double crn;
	double cin;
	double fxyr;
	double fxyi;
	double fpxyr;
	double fpxyi;
	double xp,yp;
	double den;
	double dxp;
	double dyp;
	double ds;
	double z;
	double cosz;
	double sinz;
	double rh;
	double chi;
	double dphi;
	double phi;
	long j;
	long nn;

	/* Inverse equations
	-----------------*/
	x = (x - m_falseEasting) / m_rMajor;
	y = (y - m_falseNorthing) / m_rMajor;
	xp = x;
	yp = y;
	nn = 0;

	/* Use Knuth algorithm for summing complex terms, to convert Modified-
	Stereographic conformal to Oblique Stereographic coordinates.
	--------------------------------------------------------------------*/
	do
	{
		r = xp + xp;
		s = xp * xp + yp * yp;
		ar = m_acoef[m_n];
		ai = m_bcoef[m_n];
		br = m_acoef[m_n -1];
		bi = m_bcoef[m_n - 1];
		cr = (double) (m_n) * ar;
		ci = (double) (m_n) * ai;
		dr = (double) (m_n -1) * br;
		di = (double) (m_n -1) * bi;

		for (j = 2; j <= m_n; j++)
		{
			arn = br + r * ar;
			ain = bi + r * ai;
			if (j < m_n)
			{
				br = m_acoef[m_n -j] - s * ar;
				bi = m_bcoef[m_n - j] - s * ai;
				ar = arn;
				ai = ain;
				crn = dr  + r * cr;
				cin = di  + r * ci;
				dr = (double) (m_n - j) * m_acoef[m_n -j] - s * cr;
				di = (double) (m_n - j) * m_bcoef[m_n -j] - s * ci;
				cr = crn;
				ci = cin;
			}
		}
		br = -s * ar;
		bi = -s * ai;
		ar = arn;
		ai = ain;
		fxyr = xp * ar - yp * ai + br - x;
		fxyi = yp * ar + xp * ai + bi - y;
		fpxyr = xp * cr - yp * ci + dr;
		fpxyi = yp * cr + xp * ci + ci;
		den = fpxyr * fpxyr + fpxyi * fpxyi;
		dxp = -(fxyr * fpxyr + fxyi * fpxyi) / den;
		dyp = -(fxyi * fpxyr - fxyr * fpxyi) / den;
		xp = xp + dxp;
		yp = yp + dyp;
		ds = fabs(dxp) + fabs(dyp);
		nn++;
		if (nn > 20)
		{
			setError(235);
			return;
		}
	}
	while (ds > EPSLN);

	/* convert Oblique Stereographic coordinates to LAT/LONG
	------------------------------------------------------*/
	rh = sqrt(xp * xp + yp * yp);
	z = 2.0 * atan(rh / 2.0);
	Util::gctp_sincos(z,&sinz,&cosz);
	m_longitude = m_centerLon;
	if (fabs(rh) <= EPSLN)
	{
		m_latitude = m_centerLat;
		return;
	}
	chi = Util::asinz(cosz * m_sinCenterLat + (yp * sinz * m_cosCenterLat) / rh);
	nn = 0;
	phi = chi;
	do
	{
		esphi = m_e * sin(phi);
		dphi = 2.0 * atan(tan((HALF_PI + chi) / 2.0) * 
			pow(((1.0 + esphi) / (1.0 - esphi)),(m_e / 2.0))) - HALF_PI - phi;
		phi += dphi;
		nn++;
		if (nn > 20)
		{
			setError(236);
			return;
		}
	}
	while(fabs(dphi) > EPSLN);

	m_latitude = phi;
	m_longitude = Util::adjust_lon (m_centerLon + atan2((xp * sinz), (rh * m_cosCenterLat * cosz - yp *
					m_sinCenterLat * sinz)));
}