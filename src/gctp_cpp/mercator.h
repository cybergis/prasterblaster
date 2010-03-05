#ifndef MERCATOR_H
#define MERCATOR_H
#include "projection.h"

//!This is the object used for the Mercator projection.
class Mercator: public Projection 
{
public:
	Mercator();

	Mercator(double gctpParams[],  ProjUnit units, ProjDatum dat);
	
protected:

	//!Eccentricity
	double m_e;

	//!Eccentricity squared.
	double m_es;

	double m_m1;

	//! See documentation for Projection
	void _init();

	//! See documentation for Projection
	void _forward(double lon, double lat);

	//! See documentation for Projection
	void _inverse(double x, double y);


};

#endif
	