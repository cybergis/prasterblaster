#ifndef PROJECTION_H
#define PROJECTION_H

#include <math.h>
#include <string>
#include "util.h"
#include "constants.h"

//! Base projection class.
/*! This class provides an interface that all inheriting projection objects
	must adhere to.
*/
class Projection
{
   
  public:

    Projection();

	//! Constructor 
	/*! This constructor allows the user to construct a fully initialized projection 
		object that is ready to be used.
		\param gctpParameters The 15 required projection parameters
		\param units The units that this projection will use
		\param datum The datum that this projection uses
		\param spheroid The spheroid that this projection uses
		*/
	Projection (double gctpParameters[], ProjUnit units, ProjDatum dat);
    
	//! Perform a forward transformation.

	/*! This function transforms input lat/lon coordinates to the desired coordinate system.
		If desired the user can use non-NULL pointers for x and y to obtain the output
		coordinates, but member variables are also set with these values that can be
		obtained using the x() and y() member functions.
		\param lon Input Longitude
		\param lat Input Latitude
		\param x Optional storage for output x coordinate
		\param y Optional storage for output y coordinate
	*/
	void forward ( double lon, double lat, double* x = NULL, double* y = NULL );
    
	//! Perform an inverse transformation.
	/*! This function transforms input x/y coordinates to the geographic coordinate system.
		If desired the user can use non-NULL pointers for lat and lon to obtain the output
		coordinates, but member variables are also set with these values that can be
		obtained using the lat() and lon() member functions.
		\param x Input x cooridnate
		\param y Input y coordinate
		\param lat Optional storage for output latitude
		\param lon Optional storage for output longitude
	*/
	void inverse ( double x, double y, double* lon = NULL, double* lat = NULL );

	//! Get the x coordinate from the most recent forward() call.
	double x () {return m_x_coord;}

	//! Get the y coordinate from the most recent forward() call.
	double y () {return m_y_coord;}

	//! Get the latitude from the most recent inverse() call.
	double lat () {return m_latitude;}

	//! Get the longitude from the most recent inverse() call.
	double lon () {return m_longitude;}

	//! Get both the x and y coordinates from the most recent forward() call.
    void xy (double* x, double* y);

	//! Get both the lat and lon from the most recent inverse() call.
    void latLon ( double* lat, double* lon );

	//! Get the gctp paremeter array.
	double* params() {return m_gctpParams;}

	//! Get a specific parameter in the gctp parameter array.
	double param(size_t index);

	//! Get the name of the projection.
	std::string name() {return m_name;}

	//! Get the projection number.
	ProjCode number() { return m_number;}

	//! Get the units being used.
	ProjUnit units() {return m_unitCode;}

	//! Get the spheroid being used.
	ProjDatum datum() {return m_datum;}

	//! Set the units being used.
	void setUnits(ProjUnit units) {m_unitCode = units; setInit();}
	
	//! Set the spheroid being used.
	void setDatum(ProjDatum dat) {m_datum = dat; setInit(); setRadii();}
		
	//! Set the false easting.
	void setFE(double fe) {m_falseEasting = fe; setInit();}
	
	//! Set the false northing.
	void setFN(double fn) {m_falseNorthing = fn; setInit();}
	
	//! Set sphere radii according to spheroid and parameter array
	void setRadii() {Util::sphdz(m_datum, m_gctpParams, &m_rMajor, &m_rMinor, &m_radius);
					 setInit();}
	
	//! Set the major sphere radius
	void setRMajor(double rMajor) {m_rMajor = rMajor; setInit();}
	
	//! Set the minor sphere radius
	void setRMinor(double rMinor) {m_rMinor = rMinor; setInit();}

	//! Set the radius of the sphere.
	void setRadius(double radius) {m_radius = radius; setInit();}

	//!Set the center longtitude of the projection.
	/*! This function sets the center longitude for
		the projection.
		\param lon The center longitude in decimal degrees.
	*/
	void setCenterLon(double lon);

	//!Set the center latitude of the projection.
	/*! This function sets the center latitude for
		the projection.
		\param lat The center latitude in decimal degrees.
	*/
	void setCenterLat(double lat);
	
	//!Set the latitude of the first standard parallel.
	/*! This function sets the latitude of the first 
		standard parallel.
		\param lat The latitude of the first standard parallel in decimal degrees.
	*/
	void setStdParallelLat1(double lat);

	//!Set the latitude of the second standard parallel.
	/*! This function sets the latitude of the second 
		standard parallel.
		\param lat The latitude of the second standard parallel in decimal degrees.
	*/
	void setStdParallelLat2(double lat);

	//! Get the latitude of the first standard parallel.
	double stdParallelLat1() {return m_stdParallelLat1;}

	//! Get the latitude of the second standard parallel.
	double stdParallelLat2() {return m_stdParallelLat2;}

	//! Set the gctp parameter array.
	void setParams(double gctpParams[]);

	//! Set a particular parameter in the gctp parameter array.
	void setParam(size_t index, double value);

	//! Get the current error code.
	long error() {return m_errorCode;}

	//! Check if an error has occured.
	bool errorOccured() {return(m_errorCode != 0);}


protected:
	
	//! Set the error flag to indicate that an error has occured.
	void setError(long errorCode) {m_errorCode = errorCode;}

	//! Clear the error flag.
	void clearError() {m_errorCode = 0;}
	
	//! Numeric error code of latest error
	long m_errorCode;

	//! The name of the projection.
	std::string m_name;

	//! The numeric identifier of the projection.
    ProjCode m_number;
    
	//! The numeric identifier for the units of this projection.
	ProjUnit m_unitCode;

	//! The numeric identifier for the spheroid of this projection.
    ProjDatum m_datum;

	//! The longitude value produced from an inverse transformation.
    double m_longitude;

	//! The latitude value produced from an inverse transformation.
	double m_latitude;

	//! The x coordinate produced from a forward transformation.
    double m_x_coord;

	//! The y coordinate produced from a forward transformation.
    double m_y_coord;

	//! False easting value.
	double m_falseEasting;

	//! False northing value.
	double m_falseNorthing;

	//! Major radius of sphere.
	double m_rMajor;

	//! Minor radius of sphere.
	double m_rMinor;

	//! Radius of sphere.
	double m_radius;

	//!Center longitude / longitude of central meridian of the projection.
	double m_centerLon;

	//!Center latitude / latitude of central merridian of the projection.
	double m_centerLat;

	//! Latitude of first standard parallel.
	double m_stdParallelLat1;

	//! Latitude of second standard parallel.
	double m_stdParallelLat2;

	//! Array of 15 projection parameters (as used in the original GCTP).
    double m_gctpParams[COEFCT];

	//!Do we need to do an initialization?
	bool m_initNeeded;

	//!Do we need to reload the members from the parameter array?
	bool m_paramLoadNeeded;

	//!Perform all intializations needed for forward and inverse transformations.
	/*! This function calls the virtual function _init() which is overriden
		by each derived class and performs other setup as needed by every
		projection.
	*/
	void init();

	//!This function performs the actual initialization for each projection.
	virtual void _init() = 0;

	//! Set the number of the projection.
	void setNumber(ProjCode number) {m_number = number;}

	//! Set the name of the projection.
	void setName(std::string name) {m_name = name;}

	//! Toggle forward and inverse initialization flags.
	void setInit(bool set = true) {m_initNeeded = set;}

	//! Do we need to reinitialize?
	bool initNeeded() {return m_initNeeded;}
	
	//!Toggle reload from parameter array.
	void setParamLoad(bool set = true) {m_paramLoadNeeded = set; setInit();}

	//!Do we need to reload the members from the parameter array?
	bool paramLoadNeeded() {return(m_paramLoadNeeded);}

	//!Load member variables with their corresponding values in the parameter array.
	/*! If any derived class needs to load parameters from the gctp parameter array
		aside from false easting, false northing, and sphere radii this function 
		must be overloaded.
	*/
	virtual void _loadFromParams();

	//!Converts the packed DMS angle specified by "angle" to radians and assigns it to "member".
	/*! This function was created to avoid excessive code repetition, as the majority
		of projections need to convert specific parameters from the packed DMS angle 
		format to radians and store that value in a particular member variable.
	*/
	void convertAndSetAngle(double& member, double angle);

	//! Performs the forward transformation.
	/*! This function is responsible for performing
		the projection's forward transformation.
		\param lon Input longitude in decimal degrees
		\param lat Input latitude in decimal degrees
	*/
	virtual void _forward(double lon, double lat) = 0;

	//! Performs the inverse transformation.
	/*! This function is responsible for performing
		the projection's inverse transformation.
		Note that the units of the input coordinates
		must be set either in the constructor or by
		using the setUnits() function. For a list of 
		unit codes please 
		\param x Input x coordinate
		\param y Input y coordinate
	*/
	virtual void _inverse(double x, double y) = 0;

	void loadFromParams();



};

#endif
