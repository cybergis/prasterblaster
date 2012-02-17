#include "sinusoidal.h"

#include <stdio.h>
#include <math.h>

#include <ogr_spatialref.h>

Sinusoidal::Sinusoidal(): Projection()
{
	setName("Sinusoidal");
	setNumber(SNSOID);
}

Sinusoidal::Sinusoidal( double gctpParameters[15], ProjUnit units, ProjDatum dat) 
: Projection( gctpParameters, units, dat )
{
  setName("Sinusoidal");
  setNumber(SNSOID);
}

void Sinusoidal::_forward (double lon, double lat)
{
  double deltaLon;	/* Delta longitude (Given longitude - center */
  
  /* Forward equations */
  deltaLon = Util::adjust_lon(lon - m_centerLon);

  m_x_coord = m_radius * deltaLon * cos( lat ) + m_falseEasting;
  m_y_coord = m_radius * lat + m_falseNorthing;

}

void Sinusoidal::_inverse(double x, double y)
{
  double temp;		/* Re-used temporary variable */
  
  x -= m_falseEasting;
  y -= m_falseNorthing;

  m_latitude = y / m_radius;

  if( fabs( m_latitude ) > HALF_PI) 
  {
		setError(164);
		return;
  }

  temp = fabs( m_latitude ) - HALF_PI;
  if( fabs( temp ) > EPSLN )
  {
     temp = m_centerLon + x / (m_radius * cos( m_latitude ));
	 m_longitude = Util::adjust_lon( temp );
  }
  else
  {
     m_longitude = m_centerLon;
  }

}

void Sinusoidal::_init(  )
{
	return;
}

std::string Sinusoidal::wkt()
{
		OGRSpatialReference sr;
	int epsg = DATUM2EPSG[datum()];
	char *wkt = 0;
	std::string output = "";
	OGRErr err;


	if (epsg != -1) {
		err = sr.importFromEPSG(epsg);
		if (err != OGRERR_NONE) {
			fprintf(stderr, "Error setting EPSG\n");
		}
		err = sr.SetMollweide(param(4), param(6), param(7));
		if (err != OGRERR_NONE) {
			fprintf(stderr, "Error setting projection\n");
		}
		sr.Fixup();
		err = sr.Validate();
		err = sr.exportToPrettyWkt(&wkt);

	} else {
		return output;
	}


	if (err == OGRERR_NONE) {

		output = wkt;
		OGRFree(wkt);
	} else {
		printf("WKT Broken!\n");
		if (err == OGRERR_UNSUPPORTED_SRS) {
			printf("Unsupported SRS!\n");
			
		} else if (err == OGRERR_CORRUPT_DATA) {
			printf("SRS not well formed!\n");
		}
		output = wkt;
		OGRFree(wkt);
	}

	return output;


}




