
#ifndef INTMOLLWEIDE_H
#define INTMOLLWEIDE_H

#include "projection.h"

//!This is the object used for the Interrupted Mollweide projection.
class IntMollweide : public Projection
{
public:
	IntMollweide();
	IntMollweide(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:
	
	//!Longitude of central meridian for each region. 
	double m_centerLons[6];

	//!False easting value for each region.
	double m_falseEastings[6];

	//! See documentation for Projection
	void _init();

	//! See documentation for Projection
	void _forward(double lon, double lat);

	//! See documentation for Projection
	void _inverse(double x, double y);

};

#endif
