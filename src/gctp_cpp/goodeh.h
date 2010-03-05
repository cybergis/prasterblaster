
#ifndef GOODE_H
#define GOODE_H
#include "projection.h"

//!This is the object used for the the interrupted Goode's homosline projection.
class GoodeH : public Projection
{
public:

	GoodeH();
	GoodeH(double gctpParams[], ProjUnit units, ProjDatum dat);

	//!Get the array of longitudes of central meridians for the 12 regions.
	double* centerLons() {return m_centerLons;}

	//!Get the array of false easting values for the 12 regions.
	double* falseEastings() {return m_falseEastings;}
private:

	//!Longitude of central meridians for each of the 12 regions.
	double m_centerLons[12];

	//!False easting values for each of the 12 regions.
	double m_falseEastings[12];

	//!See documentation for projection.
	void _init();

	//!See documentation for projection.
	void _forward(double lon, double lat);

	//!See documentation for projection.
	void _inverse(double x, double y);
};

#endif