#include "robinson.h"

Robinson::Robinson() : Projection()
{        

	setNumber(ROBIN);
	setName("Robinson");

	m_pr[1]= -0.062; 
	m_xlr[1]=0.9986; 
	m_pr[2]=0.0; 
	m_xlr[2]=1.0; 
	m_pr[3]=0.062; 
	m_xlr[3]=0.9986;
	m_pr[4]=0.124;   
	m_xlr[4]=0.9954;  
	m_pr[5]=0.186;
	m_xlr[5]=0.99; 
	m_pr[6]=0.248;  
	m_xlr[6]=0.9822; 
	m_pr[7]=0.31;     
	m_xlr[7]=0.973;    
	m_pr[8]=0.372;      
	m_xlr[8]=0.96;
	m_pr[9]=0.434; 
	m_xlr[9]=0.9427;
	m_pr[10]=0.4958;
	m_xlr[10]=0.9216;
	m_pr[11]=0.5571;  
	m_xlr[11]=0.8962;  
	m_pr[12]=0.6176;
	m_xlr[12]=0.8679;
	m_pr[13]=0.6769; 
	m_xlr[13]=0.835;  
	m_pr[14]=0.7346;
	m_xlr[14]=0.7986;
	m_pr[15]=0.7903;  
	m_xlr[15]=0.7597;  
	m_pr[16]=0.8435;
	m_xlr[16]=0.7186;
	m_pr[17]=0.8936; 
	m_xlr[17]=0.6732; 
	m_pr[18]=0.9394; 
	m_xlr[18]=0.6213;
	m_pr[19]=0.9761;  
	m_xlr[19]=0.5722;  
	m_pr[20]=1.0;  
	m_xlr[20]=0.5322;

	for (int i = 0; i < 21; i++)
		m_xlr[i] *= 0.9858;
}

Robinson::Robinson(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat)
{
	setNumber(ROBIN);
	setName("Robinson");

	m_pr[1]= -0.062; 
	m_xlr[1]=0.9986; 
	m_pr[2]=0.0; 
	m_xlr[2]=1.0; 
	m_pr[3]=0.062; 
	m_xlr[3]=0.9986;
	m_pr[4]=0.124;   
	m_xlr[4]=0.9954;  
	m_pr[5]=0.186;
	m_xlr[5]=0.99; 
	m_pr[6]=0.248;  
	m_xlr[6]=0.9822; 
	m_pr[7]=0.31;     
	m_xlr[7]=0.973;    
	m_pr[8]=0.372;      
	m_xlr[8]=0.96;
	m_pr[9]=0.434; 
	m_xlr[9]=0.9427;
	m_pr[10]=0.4958;
	m_xlr[10]=0.9216;
	m_pr[11]=0.5571;  
	m_xlr[11]=0.8962;  
	m_pr[12]=0.6176;
	m_xlr[12]=0.8679;
	m_pr[13]=0.6769; 
	m_xlr[13]=0.835;  
	m_pr[14]=0.7346;
	m_xlr[14]=0.7986;
	m_pr[15]=0.7903;  
	m_xlr[15]=0.7597;  
	m_pr[16]=0.8435;
	m_xlr[16]=0.7186;
	m_pr[17]=0.8936; 
	m_xlr[17]=0.6732; 
	m_pr[18]=0.9394; 
	m_xlr[18]=0.6213;
	m_pr[19]=0.9761;  
	m_xlr[19]=0.5722;  
	m_pr[20]=1.0;  
	m_xlr[20]=0.5322;

	for (int i = 0; i < 21; i++)
		m_xlr[i] *= 0.9858;
}

void Robinson::_inverse(double x, double y)
{
	double yy;
	double p2;
	double u,v,t,c;
	double phid;
	double y1;
	long ip1;
	long i;


	/* Inverse equations
	-----------------*/
	x -= m_falseEasting;
	y -= m_falseNorthing;

	yy = 2.0 * y / PI / m_radius;
	phid = yy * 90.0;
	p2 = fabs(phid / 5.0);
	ip1 = (long) (p2 - EPSLN);
	if (ip1 == 0)
	ip1 = 1;

	/* Stirling's interpolation formula as used in forward transformation is 
	reversed for first estimation of LAT. from rectangular coordinates. LAT.
	is then adjusted by iteration until use of forward series provides correct 
	value of Y within tolerance.
	---------------------------------------------------------------------------*/
	for (i = 0;;)
	{
		u = m_pr[ip1 + 3] - m_pr[ip1 + 1];
		v = m_pr[ip1 + 3] - 2.0 * m_pr[ip1 + 2] + m_pr[ip1 + 1];
		t = 2.0 * (fabs(yy) - m_pr[ip1 + 2]) / u;
		c = v / u;
		p2 = t * (1.0 - c * t * (1.0 - 2.0 * c * t));

		if ((p2 >= 0.0) || (ip1 == 1))
			{
			if (y >= 0)
				phid = (p2 + (double) ip1 ) * 5.0;
			else
				phid = -(p2 + (double) ip1 ) * 5.0;
		      
			do
			{
				p2 = fabs(phid / 5.0);
				ip1 = (long) (p2 - EPSLN);
				p2 -= (double) ip1;
		 
				if (y >= 0)
				y1 = m_radius * (m_pr[ip1 +2] + p2 *(m_pr[ip1 + 3] - m_pr[ip1 +1]) / 2.0 + p2 
						* p2 * (m_pr[ip1 + 3] - 2.0 * m_pr[ip1 + 2] + m_pr[ip1 + 1])/2.0) 
						* PI / 2.0; 
				else
				y1 = -m_radius * (m_pr[ip1 +2] + p2 *(m_pr[ip1 + 3] - m_pr[ip1 +1]) / 2.0 + p2 
						* p2 * (m_pr[ip1 + 3] - 2.0 * m_pr[ip1 + 2] + m_pr[ip1 + 1])/2.0) 
						* PI / 2.0; 
				phid += -180.0 * (y1 - y) / PI / m_radius;
				i++;
				if (i > 75)
				{
					setError(234);
					return;
				}
			}
		while (fabs(y1 - y) > .00001);
		break;
		}
	else
	{
		ip1 -= 1;
		if (ip1 < 0)
		{
			setError(234);
			return;
		}
	}
}
	m_latitude  = phid * .01745329252;

	/* calculate  LONG. using final LAT. with transposed forward Stirling's 
	interpolation formula.
	---------------------------------------------------------------------*/
	m_longitude = m_centerLon + x / m_radius / (m_xlr[ip1 + 2] + p2 * (m_xlr[ip1 + 3] - m_xlr[ip1 + 1])
						/ 2.0 + p2 * p2 * (m_xlr[ip1 + 3] - 2.0 * m_xlr[ip1 + 2] + 
						m_xlr[ip1 + 1]) / 2.0);
	m_longitude = Util::adjust_lon(m_longitude);
}

void Robinson::_forward(double lon, double lat)
{
	double dlon;
	double p2;
	long ip1;

	/* Forward equations
	-----------------*/
	dlon = Util::adjust_lon(lon - m_centerLon);
	p2 = fabs(lat / 5.0 / .01745329252);
	ip1 = (long) (p2 - EPSLN);

	/* Stirling's interpolation formula (using 2nd Diff.)
	---------------------------------------------------*/
	p2 -= (double) ip1;
	m_x_coord = m_radius * (m_xlr[ip1 + 2] + p2 * (m_xlr[ip1 + 3] - m_xlr[ip1 + 1]) / 2.0 +
			p2 * p2 * (m_xlr[ip1 + 3] - 2.0 * m_xlr[ip1 + 2] + m_xlr[ip1 + 1])/2.0) * 
			dlon + m_falseEasting;

	if (lat >= 0)
		m_y_coord = m_radius * (m_pr[ip1 + 2] + p2 * (m_pr[ip1 + 3] - m_pr[ip1 +1]) / 2.0 + p2 * p2 *
					(m_pr[ip1 + 3] - 2.0 * m_pr[ip1 + 2] + m_pr[ip1 + 1]) / 2.0) * PI / 2.0 +
				m_falseNorthing;
	else
		 m_y_coord = -m_radius * (m_pr[ip1 + 2] + p2 * (m_pr[ip1 + 3] - m_pr[ip1 +1]) / 2.0 + p2 * p2 *
             (m_pr[ip1 + 3] - 2.0 * m_pr[ip1 + 2] + m_pr[ip1 + 1]) / 2.0) * PI / 2.0 + m_falseNorthing;
}

void Robinson::_init()
{
	return;
}
