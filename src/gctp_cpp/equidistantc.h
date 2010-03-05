
#ifndef EQUID_C_H
#define EQUID_C_H
#include "projection.h"

//!This is the object used for the equidistant conic projection.
class EquidistantC : public Projection
{
public:

	EquidistantC();

	EquidistantC(double gctpParams[], ProjUnit units, ProjDatum dat);

	//!Set the mode of the projection.
	/*! This function allows you to set the mode
		of the projection.
		\param mode The mode of the projection (0 = A, 1 = B)
	*/
	void setMode(int mode) {m_mode = mode; setInit();}

	int mode() {return m_mode;}

protected:
	
	//!Eccentricity.
	double m_e;
	
	//!Eccentricity squared.
	double m_es;
	
	//!Eccentricity constant.
	double m_e0;

	//!Eccentricity constant.
	double m_e1;

	//!Eccentricity constant.	
	double m_e2;
	
	//!Eccentricity constant.
	double m_e3;
	
	//!Small value m.
	double m_ml0;
	
	double m_ns;
	
	double m_g;
	
	double m_rh;

	int m_mode;
	
	//! See documentation for Projection
	void _init();

	//! See documentation for Projection
	void _loadFromParams();

	//! See documentation for Projection
	void _forward(double lon, double lat);

	//! See documentation for Projection
	void _inverse(double x, double y);

};

#endif

