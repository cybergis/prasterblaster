#ifndef TRANS_MERC_H
#define TRANS_MERC_H
#include "projection.h"

//! This is the object used for the Transverse Mercator projection
class TransverseMercator : public Projection 
{
public:
	TransverseMercator();
	TransverseMercator(double gctpParams[], ProjUnit units, ProjDatum datum);

	//! Set the projection scale factor.
	void setScaleFactor(double scale) {m_scaleFactor = scale; setInit();}

	//! Get the scale factor.
	double scaleFactor() {return m_scaleFactor;}

protected:

	//! See documentation for Projection
	void _init();

	//! See documentation for Projection
	void _loadFromParams();

	//! See documentation for Projection
	void _forward(double lon, double lat);

	//! See documentation for Projection
	void _inverse(double x, double y);

	double m_scaleFactor;

	double m_e0;
	double m_e1;
	double m_e2;
	double m_e3;
	double m_e;
	double m_es;
	double m_esp;
	double m_ml0;
	double m_ind;
};

#endif

	