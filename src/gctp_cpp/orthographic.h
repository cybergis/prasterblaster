
#ifndef ORTHO_H
#define ORTHO_H

#include "projection.h"

//! This is the object used for the Orthographic projection.
class Orthographic : public Projection
{
public:
	Orthographic();
	Orthographic(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	//! Sin of the senter latitude of the projection.
	double m_sinCenterLat;

	//! Cosine of the center latitude of the projection.
	double m_cosCenterLat;

	//! See documentation for Projection
	void _init();

	//! See documentation for Projection
	void _inverse(double x, double y);

	//! See documentation for Projection
	void _forward(double lon, double lat);

};

#endif
