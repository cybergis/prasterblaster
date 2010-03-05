
#ifndef SPACE_OB_MERC_H
#define SPACE_OB_MERC_H

#include "projection.h"

//! This is the object used for the Space Oblique Mercator projection.
class SpaceObMerc : public Projection
{
public:

	SpaceObMerc();
	SpaceObMerc(double gctpParams[], ProjUnit units, ProjDatum dat);

	//! Set the satellite number
	void setSatNum(long num) {m_satNum = num; setInit();}

	//! Set the landsat path number
	void setPath(long path) {m_path = path; setInit();}

	//! Set the mode (0=A, 1=B)
	void setMode(long mode) {m_mode = mode; setInit();}

	//! Set the inclination of orbit at the ascending node, counter-clockwise from the equator.
	void setAlf(double val);
	
	//! Set the period of the satellite revolution, in minutes.
	void setTime(double time) {m_time = time; setInit();}

	//! Set the end of path flag, (0 = start, 1 = end).
	void setStart(double start) {m_start = start; setInit();}

	//! Get the landsat satellite number.
	long satNum() {return m_satNum;}

	//! Get the landsat path number.
	long path() {return m_path;}

	//! Get the mode of the projection, (0=A, 1=B).
	long mode() {return m_mode;}

	//! Set the inclination of orbit at the ascending node, counter-clockwise from the equator.
	double alf() {return m_alf;}

	//! Set the period of the satellite revolution, in minutes.
	double time() {return m_time;}

	//! Set the end of path flag, (0 = start, 1 = end).
	double start() {return m_start;}

protected:

	long m_satNum;
	long m_path;
	long m_mode;
	double m_alf;
	double m_time;
	double m_a;
	double m_b;
	double m_a2;
	double m_a4;
	double m_c1;
	double m_c3;
	double m_q;
	double m_t;
	double m_u;
	double m_w;
	double m_xj;
	double m_p21;
	double m_sa;
	double m_ca;
	double m_es;
	double m_s;
	double m_start;
	double m_centerLon;

	//! Calculate a,b, and c coefficients to convert from transform lat/lon to SOM rectangular coordinates.
	void som_series(double* fb, double* fa2, double* fa4, double* fc1, double* fc3, double* dlam);

	//! See documentation for Projection.
	void _init();

	//! See documentation for Projection.
	void _loadFromParams();

	//! See documentation for Projection.
	void _inverse(double x, double y);

	//! See documentation for Projection.
	void _forward(double lon, double lat);
};

#endif