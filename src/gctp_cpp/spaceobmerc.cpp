#include "spaceobmerc.h"

SpaceObMerc::SpaceObMerc() : Projection(),
m_satNum(0), m_path(0), m_mode(0),
m_alf(0.0), m_time(0.0),
m_a(0.0), m_b(0.0),
m_a2(0.0), m_a4(0.0), m_c1(0.0),
m_c3(0.0), m_q(0.0), m_t(0.0),
m_u(0.0), m_w(0.0), m_xj(0.0),
m_p21(0.0), m_sa(0.0), m_ca(0.0),
m_es(0.0), m_s(0.0), m_start(0.0),
m_centerLon(0.0)

{
	setNumber(SOM);
	setName("Space Oblique Mercator");
}

SpaceObMerc::SpaceObMerc(double gctpParams[], ProjUnit units, ProjDatum dat):
Projection(gctpParams, units, dat),
m_satNum(0), m_path(0), m_mode(0),
m_alf(0.0), m_time(0.0),
m_a2(0.0), m_a4(0.0), m_c1(0.0),
m_c3(0.0), m_q(0.0), m_t(0.0),
m_u(0.0), m_w(0.0), m_xj(0.0),
m_p21(0.0), m_sa(0.0), m_ca(0.0),
m_es(0.0), m_s(0.0), m_start(0.0),
m_centerLon(0.0)

{
	setNumber(SOM);
	setName("Space Oblique Mercator");
	setParamLoad();
}

void SpaceObMerc::_init()
{
	long i;
	double alf,e2c,e2s,one_es;
	double dlam,fb,fa2,fa4,fc1,fc3,suma2,suma4,sumc1,sumc3,sumb;

	m_a = m_rMajor;
	m_b = m_rMinor;
	m_es = 1.0 - SQUARE(m_rMinor/m_rMajor);

	if (m_mode != 0)
	{
		alf = m_alf;
		m_p21 = m_time/1440.0;
	}
	else
	{
		if (m_satNum < 4)
		{
			alf = 99.092 * D2R;
			m_p21=103.2669323/1440.0;
			m_centerLon = (128.87 - (360.0/251.0 * m_path)) * D2R;
		}
		else
		{
			alf = 98.2 * D2R;
			m_p21 = 98.8841202/1440.0;
			m_centerLon = (129.30 - (360.0/233.0 * m_path)) * D2R;
		}
		m_start=0.0;
	}


	m_ca=cos(alf);

	if (fabs(m_ca)<1.e-9) 
		m_ca=1.e-9;

	m_sa=sin(alf);
	e2c = m_es*m_ca*m_ca;
	e2s = m_es*m_sa*m_sa;
	m_w = (1.0-e2c)/(1.0-m_es);
	m_w = m_w*m_w-1.0;
	one_es = 1.0-m_es;
	m_q = e2s / one_es;
	m_t = (e2s*(2.0-m_es)) / (one_es*one_es);
	m_u = e2c / one_es;
	m_xj = one_es*one_es*one_es;
	dlam = 0.0;
	som_series(&fb,&fa2,&fa4,&fc1,&fc3,&dlam);
	suma2 = fa2;
	suma4 = fa4;
	sumb = fb;
	sumc1 = fc1;
	sumc3 = fc3;
	for(i=9;i<=81;i+=18)
	{
		dlam = i;
		som_series(&fb,&fa2,&fa4,&fc1,&fc3,&dlam);
		suma2 = suma2+4.0*fa2;
		suma4 = suma4+4.0*fa4;
		sumb = sumb+4.0*fb;
		sumc1 = sumc1+4.0*fc1;
		sumc3 = sumc3+4.0*fc3;
	}
	for(i=18; i<=72; i+=18)
	{
		dlam = i;
		som_series(&fb,&fa2,&fa4,&fc1,&fc3,&dlam);
		suma2 = suma2+2.0*fa2;
		suma4 = suma4+2.0*fa4;
		sumb = sumb+2.0*fb;
		sumc1 = sumc1+2.0*fc1;
		sumc3 = sumc3+2.0*fc3;
	}

	dlam = 90.0;
	som_series(&fb,&fa2,&fa4,&fc1,&fc3,&dlam);
	suma2 = suma2+fa2;
	suma4 = suma4+fa4;
	sumb = sumb+fb;
	sumc1 = sumc1+fc1;
	sumc3 = sumc3+fc3;
	m_a2 = suma2/30.0;
	m_a4 = suma4/60.0;
	m_b = sumb/30.0;
	m_c1 = sumc1/15.0;
	m_c3 = sumc3/45.0;


}

void SpaceObMerc::som_series(double* fb, double* fa2, double* fa4, double* fc1, double* fc3, double* dlam)
{

	double sd,sdsq,h,sq,fc;

	*dlam= *dlam*0.0174532925;               /* Convert dlam to radians */
	sd=sin(*dlam);
	sdsq=sd*sd;
	m_s=m_p21*m_sa*cos(*dlam)*sqrt((1.0+m_t*sdsq)/((1.0+m_w*sdsq)*(1.0+m_q*sdsq)));
	h = sqrt((1.0+m_q*sdsq)/(1.0+m_w*sdsq))*(((1.0+m_w*sdsq)/((1.0+m_q*sdsq)*(1.0+
		m_q*sdsq)))-m_p21*m_ca);
	sq=sqrt(m_xj*m_xj+m_s*m_s);
	*fb=(h*m_xj-m_s*m_s)/sq;
	*fa2= *fb*cos(2.0* *dlam);
	*fa4= *fb*cos(4.0* *dlam);
	fc=m_s*(h+m_xj)/sq;
	*fc1=fc*cos(*dlam);
	*fc3=fc*cos(3.0* *dlam);
}

void SpaceObMerc::_loadFromParams()
{
	Projection::_loadFromParams();
	setPath((long)m_gctpParams[3]);
	setSatNum((long)m_gctpParams[2]);
	if(m_gctpParams[12] == 0) 
	{
		setMode(1);
		setAlf(m_gctpParams[3]);
		setTime(m_gctpParams[8]);
		setStart(m_gctpParams[10]);
	}

	else
		setMode(0);

}

void SpaceObMerc::setAlf(double val)
{
	long err = 0;
	err = Util::DMSToRad(val);
	if(err != 0)
	{
		setError(err);
		return;
	}
	m_alf = val;

	setInit();
}


void SpaceObMerc::_inverse(double x, double y)
{
	double tlon,conv,sav,sd,sdsq,blon,dif,st,defac,actan,tlat,dd,bigk,bigk2,xlamt;
	double sl,scl,dlat,dlon;
	long inumb;

	/* Inverse equations. Begin inverse computation with approximation for tlon. 
	Solve for transformed long.
	---------------------------*/
	x -= m_falseEasting; 
	y -= m_falseNorthing;

	tlon= x/(m_a*m_b);
	conv=1.e-9;
	for(inumb=0;inumb<50;inumb++)
	{
		sav=tlon;
		sd=sin(tlon);
		sdsq=sd*sd;
		m_s=m_p21*m_sa*cos(tlon)*sqrt((1.0+m_t*sdsq)/((1.0+m_w*sdsq)*(1.0+m_q*sdsq)));
		blon=(x/m_a)+(y/m_a)*m_s/m_xj-m_a2*sin(2.0*tlon)-m_a4*sin(4.0*tlon)-(m_s/m_xj)*
			(m_c1*sin(tlon)+m_c3*sin(3.0*tlon)); 
		tlon=blon/m_b;
		dif=tlon-sav;
		if(fabs(dif)<conv)
			break; 
	}
	if(inumb>=50)  
	{
		setError(214);
		return;
	}

	/* Compute transformed lat.
	------------------------*/
	st=sin(tlon);
	defac=exp(sqrt(1.0+m_s*m_s/m_xj/m_xj)*(y/m_a-m_c1*st-m_c3*sin(3.0*tlon)));
	actan=atan(defac);
	tlat=2.0*(actan-(PI/4.0));

	/* Compute geodetic longitude
	--------------------------*/
	dd=st*st;
	if(fabs(cos(tlon))<1.e-7) tlon=tlon-1.e-7;
	bigk=sin(tlat); 
	bigk2=bigk*bigk;
	xlamt=atan(((1.0-bigk2/(1.0-m_es))*tan(tlon)*m_ca-bigk*m_sa*sqrt((1.0+m_q*dd)
		*(1.0-bigk2)-bigk2*m_u)/cos(tlon))/(1.0-bigk2*(1.0+m_u)));

	/* Correct inverse quadrant
	------------------------*/
	if(xlamt>=0.0) sl=1.0;
	if(xlamt<0.0) sl= -1.0;
	if(cos(tlon)>=0.0) scl=1.0;
	if(cos(tlon)<0.0) scl= -1.0;
	xlamt=xlamt-((PI/2.0)*(1.0-scl)*sl);
	dlon=xlamt-m_p21*tlon;

	/* Compute geodetic latitude
	-------------------------*/
	if(fabs(m_sa)<1.e-7)dlat=asin(bigk/sqrt((1.0-m_es)*(1.0-m_es)+m_es*bigk2));
	if(fabs(m_sa)>=1.e-7)dlat=atan((tan(tlon)*cos(xlamt)-m_ca*sin(xlamt))/((1.0-m_es)*m_sa));
	m_longitude = Util::adjust_lon(dlon+m_centerLon);
	m_latitude = dlat;
}

void SpaceObMerc::_forward(double lon, double lat)
{
	long n,l;
	double delta_lon;
	double rlm,tabs,tlam,xlam,c,xlamt,ab2,ab1,xlamp,sav;
	double d,sdsq,sd,tanlg,xtan,tphi,dp,rlm2;
	double scl,tlamp,conv,delta_lat,radlt,radln;
	bool end = false;
	/* Forward equations
	-----------------*/
	conv=1.e-7;
	delta_lat = lat;
	delta_lon= lon-m_centerLon;

	/* Test for latitude and longitude approaching 90 degrees
	----------------------------------------------------*/
	if (delta_lat>1.570796) 
		delta_lat=1.570796;
	if (delta_lat<-1.570796) 
		delta_lat= -1.570796;

	radlt=delta_lat;
	radln=delta_lon;

	if(delta_lat>=0.0)
		tlamp=PI/2.0; 
	if(m_start != 0.0)
		tlamp=2.5*PI;
	if(delta_lat<0.0) 
		tlamp=1.5*PI;
	n=0;

	while(1)
	{
		sav=tlamp;
		l=0;
		end = false;
		xlamp=radln+m_p21*tlamp;
		ab1=cos(xlamp);
		
		if(fabs(ab1)<conv) 
			xlamp=xlamp-1.e-7;
		
		if(ab1>=0.0) 
			scl=1.0;
		
		if(ab1<0.0) 
			scl= -1.0;

		ab2=tlamp-(scl)*sin(tlamp)*HALF_PI;

		while(1)
		{
			xlamt=radln+m_p21*sav;
			c=cos(xlamt);

			if (fabs(c)<1.e-7) 
				xlamt=xlamt-1.e-7;

			xlam=(((1.0-m_es)*tan(radlt)*m_sa)+sin(xlamt)*m_ca)/c;
			tlam=atan(xlam);
			tlam=tlam+ab2;
			tabs=fabs(sav)-fabs(tlam);

			if(fabs(tabs)<conv) 
			{
				rlm=PI*LANDSAT_RATIO;
				rlm2=rlm+2.0*PI;
				n++;

				if(n >= 3)
					end = true;

				if((tlam>rlm) && (tlam<rlm2)) 
					end = true;

				if(tlam<rlm)
					tlamp=2.50*PI;
				if(tlam>=rlm2) 
					tlamp=HALF_PI;

				break;
			}

			l++;

			if (l > 50) 
			{
				setError(216);
				return;
			}

			sav=tlam;
		}

		if(end)
			break;
	}

	
	
	
	dp=sin(radlt);
	tphi=asin(((1.0-m_es)*m_ca*dp-m_sa*cos(radlt)*sin(xlamt))/sqrt(1.0-m_es*dp*dp));

	/* compute x and y
	---------------*/
	xtan = (PI/4.0) + (tphi/2.0);
	tanlg = log(tan(xtan));
	sd=sin(tlam);
	sdsq=sd*sd;
	m_s=m_p21*m_sa*cos(tlam)*sqrt((1.0+m_t*sdsq)/((1.0+m_w*sdsq)*(1.0+m_q*sdsq)));
	d=sqrt(m_xj*m_xj+m_s*m_s);
	m_x_coord = m_b * tlam + m_a2 * sin(2.0*tlam) + m_a4*sin(4.0*tlam)-tanlg*m_s/d;
	m_x_coord *= m_a;
	m_y_coord = m_c1*sd+m_c3*sin(3.0*tlam)+tanlg*m_xj/d;
	m_y_coord *= m_a;

	m_x_coord += m_falseEasting;
	m_y_coord += m_falseNorthing;


}

