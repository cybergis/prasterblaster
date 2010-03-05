
#ifndef MILLER_H
#define MILLER_H

#include "projection.h"

//!This is the object used for the Miller projection.
class Miller : public Projection
{
public:

	Miller();
	Miller(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	//! See documentation for Projection
	void _init();

	//! See documentation for Projection
	void _forward(double lon, double lat);

	//! See documentation for Projection
	void _inverse(double x, double y);

};

#endif