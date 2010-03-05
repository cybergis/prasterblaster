
#ifndef WAGNERVII_H
#define WAGNERVII_H

#include "projection.h"

//! This is the object used for the Wagner VII projection.
class WagnerVII : public Projection
{
public:
	WagnerVII();
	WagnerVII(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	//! See documentation for projection.
	void _init();

	//! See documentation for projection.
	void _inverse(double x, double y);

	//! See documentation for projection.
	void _forward(double lon, double lat);
};

#endif
