
#ifndef POLYCON_H
#define POLYCON_H
#include "projection.h"

//! This is the object used for the Polyconic projection.
class Polyconic: public Projection 
{
public:

	Polyconic();
	
	Polyconic(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	/* eccentricity constansts */
	double m_e;
	double m_es;
	double m_e0;
	double m_e1;
	double m_e2;
	double m_e3;
	double m_ml0; //small value m;

	//! See documentation for Projection
	void _init();

	//! See documentation for Projection
	void _forward(double lon, double lat);

	//! See documentation for Projection
	void _inverse(double x, double y);

};

#endif