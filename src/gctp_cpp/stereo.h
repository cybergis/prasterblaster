
#ifndef STEREO_H
#define STEREO_H

#include "projection.h"

//! This is the object used for the Stereographic projection.
class Stereo : public Projection
{
public:

	Stereo();
	Stereo(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	double m_sinCenterLat;
	double m_cosCenterLat;

	//! See documentation for Projection.
	void _init();

	//! See documentation for Projection.
	void _inverse(double x, double y);

	//! See documentation for Projection.
	void _forward(double lon, double lat);
};

#endif
