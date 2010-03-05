
#ifndef ROBINSON_H
#define ROBINSON_H

#include "projection.h"

//! This is the object used for the Robinson projection.
class Robinson : public Projection
{
public:

	Robinson();
	Robinson(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:
	
	double m_pr[21];
	double m_xlr[21];

	//! See documentation for Projection.
	void _init();

	//! See documentation for Projection.
	void _inverse(double x, double y);

	//! See documentation for Projection.
	void _forward(double lon, double lat);

};

#endif
