
#ifndef OBLATED_EQ_H
#define OBLATED_EQ_H

#include "projection.h"

//!This is the object used for the Oblated Equal Area Projection.
class OblatedEqArea : public Projection
{
public:

	OblatedEqArea();
	OblatedEqArea(double gctpParams[], ProjUnit units, ProjDatum dat);

	//!Set oval shape m.
	void setShapeM(double m) {m_m = m; setInit();}

	//!Set obal shape n
	void setShapeN(double n) {m_n = n; setInit();}

	//!Set oval shape angle.
	void setAngle(double theta);

	double shapeM() {return m_m;}

	double shapeN() {return m_n;}

	double angle() {return m_theta;}

protected:

	double m_theta;
	double m_m;
	double m_n;
	double m_sinLatO;
	double m_cosLatO;

	void _loadFromParams();

	void _init();

	void _forward(double lon, double lat);
	void _inverse(double x, double y);

};

#endif