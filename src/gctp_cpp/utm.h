
#ifndef UTM_H
#define UTM_H

#include "projection.h"

//! This is the object used for the UTM projection.
class UTM : public Projection
{
public:

	UTM();

	UTM(double gctpParams[], ProjUnit units, ProjDatum dat, int zone = 1);

	//! Set the projection zone.
	void setZone(int zone) {m_zone = zone;}

	//! Set the scale factor of the projection.
	void setScaleFactor(double factor) {m_scaleFactor = factor;}

	//! Get the zone of the projection.
	int zone() {return m_zone;}

	//! Get the scale factor of the projection.
	double scaleFactor() {return m_scaleFactor;}

protected:

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
	int m_zone;

	//! See documentation for projection.
	void _loadFromParams();

	//! See documentation for projection.
	void _init();

	//! See documentation for projection.
	void _forward(double lon, double lat);

	//! See documentation for projection.
	void _inverse(double x, double y);
};

#endif
