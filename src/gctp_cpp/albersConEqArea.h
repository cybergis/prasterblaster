#ifndef ALBER_CON_EQ_AREA_H
#define ALBER_CON_EQ_AREA_H
#include "projection.h"

//! This is the object used for the Albers Conical Equal-Area projection
class AlbersConEqArea: public Projection {
public:
	//! Default constructor
	AlbersConEqArea();

	//! Initialize all fields neccessary to perform our forward and inverse transformations.
	AlbersConEqArea(double gctpParameters[], ProjUnit units, ProjDatum dat);
	
protected:
  
  //! Constant C.
  double m_c;

  //! eccentricity.
  double m_e;

  //! eccentricity squared.
  double m_es;

  //! height above elipsoid.
  double m_rh;

  //! Ratio between meridians.
  double m_ns0;

  //! See documentation for Projection
  void _init();

  //! See documentation for Projection
  void _forward(double lon, double lat);

  //! See documentation for Projection
  void _inverse(double x, double y);

};

#endif