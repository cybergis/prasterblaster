#ifndef HOTIN_OB_MERC_H
#define HOTIN_OB_MERC_H

#include "projection.h"

//!This object is used for the Hotine Oblique Mercator projection.
class HotineObMerc : public Projection
{
public:

	HotineObMerc();
	HotineObMerc(double gctpParams[], ProjUnit units, ProjDatum dat);

	//!Set the scale factor for the projection.
	void setScaleFactor(double fac) {m_scaleFactor = fac; setInit();}

	//!Set the azimuthal angle for the projection.
	void setAzimuth(double angle);

	//!Set the longitude of the first point on the center line.
	void setLon1(double lon);

	//!Set the longitude of the second point on the center line.
	void setLon2(double lon);

	//!Set the latitude of the first point on the center line.
	void setLat1(double lat);

	//!Set the latitude of the second point on the center line.
	void setLat2(double lat);

	//!Set the mode of the projection.
	/*! This function allows you to set which 
		mode you will be using.
		\param mode The mode of the projection (0 = A, 1 = B)
	*/
	void setMode(int mode) {m_mode = mode; setInit();}

	//!Get the scale factor.
	double scaleFactor() {return m_scaleFactor;}

	//!Get the azimuthal angle.
	double azimuth() {return m_azimuth;}

	//!Get the longitude of the first point on the center line.
	double lon1() {return m_lon1;}

	//!Get the longitude of the second point on the center line.
	double lon2() {return m_lon2;}

	//!Get the latitude of the first point on the center line.
	double lat1() {return m_lat1;}

	//!Get the latitude of the second point on the center line.
	double lat2() {return m_lat2;}

	//!Get the mode of the projection.
	int mode() {return m_mode;}

protected:

	double m_azimuth;
	double m_scaleFactor;
	double m_e;
	double m_es;
	double m_sinCenterLat;
	double m_cosCenterLat;
	double m_bl;
	double m_al;
	double m_ts;
	double m_d;
	double m_el;
	double m_u;
	double m_singam;
	double m_cosgam;
	double m_sinaz;
	double m_cosaz;
	double m_lat1;
	double m_lat2;
	double m_lon1;
	double m_lon2;
	int m_mode;

	//!See documentation for projection.
	void _init();

	//!See documentation for projection.
	void _loadFromParams();

	//!See documentation for projection.
	void _inverse(double x, double y);
	
	//!See documentation for projection.
	void _forward(double lon, double lat);
};

#endif
