#include <stdio.h>
#include <string>

#include <ogr_spatialref.h>

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


std::string Equirectangular::wkt()
{
	std::string output = "";
	char* wkt = 0;
	OGRSpatialReference sr;
	int epsg = DATUM2EPSG[datum()];
	OGRErr err;

	if (epsg != -1) {
		err = sr.importFromEPSG(epsg);
		if (err != OGRERR_NONE) {
			fprintf(stderr, "Error setting EPSG\n");
		}
		err = sr.SetEquirectangular(0.0, param(4), param(6), param(7));
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
		printf("WKT GCTP: %s DATUM: %d\n", wkt, datum());
		output = wkt;
		OGRFree(wkt);
	}
	
	
	
	return output;
}





