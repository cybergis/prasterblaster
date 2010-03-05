
#ifndef GEN_VERT_H
#define GEN_VERT_H

#include "projection.h"

//!This is the object used for the general vertical near-side perspective projection.
class GenVertNSP : public Projection
{
public:

	GenVertNSP();

	GenVertNSP(double gctpParams[], ProjUnit units, ProjDatum dat);

	//!Set the height of the projection.
	void setHeight(double height) {m_height = height; setInit();}

	//!Get the height of the projection.
	double height() {return m_height;}

protected:

	//!Height above sphere.
	double m_p;

	//!Height of perspective point
	double m_height;

	//!Sin of the center latitude.
	double m_sinCenterLat;

	//!Cosine of the center latitude.
	double m_cosCenterLat;

	//!See documentation for projection.
	void _forward(double lon, double lat);

	//!See documentation for projection.
	void _inverse(double x, double y);

	//!See documentation for projection.
	void _init();

	//!See documentation for projection.
	void _loadFromParams();
};

#endif

