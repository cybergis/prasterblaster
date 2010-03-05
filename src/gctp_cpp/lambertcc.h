
#ifndef LAMBERT_CC_H
#define LAMBERT_CC_H
#include "projection.h"

//!This is the object used for the Lambert Conformal Conic Projection.
class LambertCC : public Projection 
{
public:
	LambertCC();
	LambertCC(double gctpParams[], ProjUnit units, ProjDatum dat);

protected:

	//! See documentation for Projection
	void _init();

	//! See documentation for Projection
	void _forward(double lon, double lat);

	//! See documentation for Projection
	void _inverse(double x, double y);

	double m_es;
	double m_e;
	double m_ns;
	double m_f0;
	double m_rh;
};
#endif
