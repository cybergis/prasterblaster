#ifndef WAGNERIV_H
#define WAGNERIV_H

#include "projection.h"

//! This is the object used for the Wagner IV projection
class WagnerIV : public Projection
{
public:

	WagnerIV();
	WagnerIV(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	//! See documentation for projection.
	void _init();

	//! See documentation for projection.
	void _inverse(double x, double y);

	//! See documentation for projection.
	void _forward(double lon, double lat);

};

#endif


