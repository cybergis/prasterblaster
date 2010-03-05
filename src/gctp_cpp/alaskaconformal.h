
#ifndef ALASKA_CON_H
#define ALASKA_CON_H

#include "projection.h"

//! This is the object used for the Alaska conformal projection.
class AlaskaConformal : public Projection
{
public:

	AlaskaConformal();
	
	AlaskaConformal(double gctpParams[], ProjUnit units, ProjDatum dat);
	
protected:
	
	//!See documentation for Projection.
	void _init();

	//!See documentation for Projection.
	void _inverse(double x, double y);

	//!See documentation for Projection.
	void _forward(double lon, double lat);
	
	double m_acoef[7];
	double m_bcoef[7];
	double m_sinCenterLat;
	double m_cosCenterLat;
	double m_e;
	int m_n;
};

#endif
