#include "projection.h"

Projection::Projection(): 
m_errorCode(0),m_longitude(0.0), m_latitude(0.0), m_x_coord(0.0), m_y_coord(0.0), 
m_falseEasting(0.0),m_falseNorthing(0.0), m_rMajor(0.0), m_rMinor(0.0), m_radius(0.0),
m_centerLon(0.0), m_centerLat(0.0), m_stdParallelLat1(0.0), m_stdParallelLat2(0.0),
m_initNeeded(false), m_paramLoadNeeded(false)

{
	for(int i = 0; i < COEFCT; i++)
		m_gctpParams[i] = 0.0;
}

Projection::Projection ( double gctpParameters[], ProjUnit units, ProjDatum dat): 
m_errorCode(0), m_unitCode(units), m_datum(dat),
m_longitude(0.0), m_latitude(0.0), m_x_coord(0.0), m_y_coord(0.0), m_falseEasting(0.0),
m_falseNorthing(0.0), m_rMajor(0.0), m_rMinor(0.0), m_radius(0.0),
m_centerLon(0.0), m_centerLat(0.0), m_stdParallelLat1(0.0), m_stdParallelLat2(0.0),
m_initNeeded(false), m_paramLoadNeeded(false)
{

  for( int index = 0; index < COEFCT; index++ )
     m_gctpParams[index] = gctpParameters[index];
  
  setParamLoad();
  return; 
}

void Projection::loadFromParams() {
	_loadFromParams();
	setParamLoad(false);
}

void Projection::_loadFromParams()
{
	setFE(m_gctpParams[6]);
	setFN(m_gctpParams[7]);
	setCenterLon(m_gctpParams[4]);
	setCenterLat(m_gctpParams[5]);
	setStdParallelLat1(m_gctpParams[2]);
	setStdParallelLat2(m_gctpParams[3]);
	setRadii();
}
void Projection::xy ( double* x, double* y )
{

  if( x != NULL )
  {
     *x = m_x_coord;
  }
  
  if( y != NULL )
  {
     *y = m_y_coord;
  }
  
}

void Projection::latLon ( double* lat, double* lon )
{
  if( lat != NULL )
  {
     *lat = m_latitude;
  }
  
  if( lon != NULL )
  {
     *lon = m_longitude;
  }
  
}

void Projection::setParams(double gctpParams[]) {
	for(int i = 0; i < COEFCT; i++) {
		m_gctpParams[i] = gctpParams[i];
	}

	setParamLoad();
	setRadii();
}

void Projection::setParam(size_t index, double value) {
	if(index < COEFCT) {
		m_gctpParams[index] = value;
		setParamLoad();
	    setRadii();
	}
}

void Projection::setCenterLat(double lat) 
{
	convertAndSetAngle(m_centerLat, lat);
	setInit();
}

void Projection::setStdParallelLat1(double lat) 
{
	convertAndSetAngle(m_stdParallelLat1, lat);
	setInit();
}

void Projection::setStdParallelLat2(double lat) 
{
	convertAndSetAngle(m_stdParallelLat2, lat);
	setInit();
}

void Projection::setCenterLon(double lon) 
{
	convertAndSetAngle(m_centerLon, lon);
	setInit();
}

double Projection::param(size_t index) {
	if(index < COEFCT)
		return m_gctpParams[index];
	else
		return(-1);
	
}

void Projection::convertAndSetAngle(double& member, double angle)
{
	double temp = angle;
	long err = 0;
	err = Util::DMSToRad(temp);
	
	if(err != 0) 
	{
		setError(err);
		return;
	}

	member = temp;
}
	

void Projection::forward(double lon, double lat, double* x, double* y)
{
	clearError();

	if(paramLoadNeeded())
		loadFromParams();

	if(errorOccured())
		return;

	if(initNeeded())
		init();

	if(errorOccured())
		return;
	
	//check if convertCoords returned an error
	setError(Util::convertCoords(DEGREE, RADIAN, lon, lat));
	if(errorOccured())
		return;

	_forward(lon, lat);

	setError(Util::convertCoords(METER, m_unitCode, m_x_coord, m_y_coord));
	if(errorOccured())
		return;
	
	if(x)
		*x = m_x_coord;
	if(y)
		*y = m_y_coord;
}

void Projection::inverse(double x, double y, double* lon, double* lat)
{
	clearError();
	
	if(paramLoadNeeded())
		loadFromParams();

	if(errorOccured())
		return;

	if(initNeeded())
		init();

	if(errorOccured())
		return;

	//check if convertCoords returned an error
	setError(Util::convertCoords(m_unitCode, METER, x, y));
	if(errorOccured())
		return;

	_inverse(x, y);

	setError(Util::convertCoords(RADIAN, DEGREE, m_longitude, m_latitude));
	
	if(errorOccured())
		return;

	if(lon)
		*lon = m_longitude;
	
	if(lat)
		*lat = m_latitude;
}

void Projection::init()
{
	_init();
	m_initNeeded = false;
}








