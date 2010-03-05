#include <stdio.h>

#include "equirectangular.h"

Equirectangular::Equirectangular() : Projection()
{
	setName("Equirectangular");
	setNumber(EQRECT);
}

Equirectangular::Equirectangular( double gctpParameters[15], ProjUnit units, ProjDatum dat) 
: Projection(gctpParameters, units, dat)
{
	setName("Equirectangular");
	setNumber(EQRECT);
}

void Equirectangular::_forward( double lon, double lat)
{
  double deltaLon;		/* delta longitude value			*/

  /* Forward equations */
  deltaLon = Util::adjust_lon( lon - m_centerLon );
  m_x_coord = m_falseEasting + m_radius * deltaLon * cos( m_centerLat );
  m_y_coord = m_falseNorthing + m_radius * lat;

}

void Equirectangular::_inverse (double x, double y)
{
  
  /* Inverse equations */
  x -= m_falseEasting;
  y -= m_falseNorthing;

  m_latitude = y / m_radius;

  if( fabs( m_latitude ) > HALF_PI ) {
	setError(174);
	return;
  }

  m_longitude = Util::adjust_lon( m_centerLon + x / ( m_radius * cos( m_centerLon )));

	 
}

void Equirectangular::_init()
{
	return;
}






