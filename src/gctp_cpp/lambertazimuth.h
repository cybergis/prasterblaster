#ifndef LAMBERT_AZ_H
#define LAMBERT_AZ_H

#include "projection.h"

//!This is the object used for the Lambert Azimuthal projection.
class LambertAzimuthal : public Projection
{
public:

	LambertAzimuthal();
	LambertAzimuthal(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	//!Sin of the center latitude.
	double m_sinCenterLat;

	//!Cosine of the center latitude.
	double m_cosCenterLat;

	//!See documentation for Projection.
	void _init();

	//!See documentation for Projection.
	void _forward(double lon, double lat);

	//!See documentation for Projection.
	void _inverse(double x, double y);

};

#endif
