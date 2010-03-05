
#ifndef GNOMONIC_H
#define GNOMONIC_H
#include "projection.h"

//!This is the object used for the gnomonic projection (***Currently Not Working!!***)
class Gnomonic : public Projection
{
public:

	Gnomonic();
	Gnomonic(double gctpParams[], ProjUnit units, ProjDatum dat);

	double sinCenterLat() {return m_sinCenterLat;}

	double cosCenterLat() {return m_cosCenterLat;}

protected:

	//!Sin of the center latitude.
	double m_sinCenterLat;

	//!Cosine of the center latitude.
	double m_cosCenterLat;

	//!See documentation for projection.
	void _init();

	//!See documentation for projection.
	void _forward(double lon, double lat);

	//!See documentation for projection.
	void _inverse(double x, double y);
};

#endif