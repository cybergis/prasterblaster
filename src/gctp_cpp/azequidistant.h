
#ifndef AZ_EQ_H
#define AZ_EQ_H

#include "projection.h"

//! This is object used for the Azimuthal Equidistant projection.
class AzEquidistant : public Projection
{
public:

	AzEquidistant();
	AzEquidistant(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	double m_sinCenterLat;
	
	double m_cosCenterLat;
	
	//! See documentation for projection.
	void _init();

	//! See documentation for projection.
	void _forward(double lon, double lat);

	//! See documentation for projection.
	void _inverse(double x, double y);

};

#endif