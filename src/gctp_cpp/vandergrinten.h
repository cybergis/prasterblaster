#ifndef VAND_H
#define VAND_H

#include "projection.h"

//! This is the object used for the Van Der Grintent projection.
class VanDerGrinten : public Projection
{
public:
	VanDerGrinten();
	VanDerGrinten(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	//! See documentation for projection.
	void _init();

	//! See documentation for projection.
	void _inverse(double x, double y);

	//! See documentation for projection.
	void _forward(double lon, double lat);
};
#endif
