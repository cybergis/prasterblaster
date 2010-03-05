

#include <boost/python.hpp>
#include <vector>

#include "constants.h"
#include "projection.h"
#include "transformer.h"

using namespace boost::python;

// Wrappers for virtual polymorphism
struct ProjectionWrap : Projection, wrapper<Projection>
{
	ProjectionWrap()
		{
			return;
		}
	ProjectionWrap(boost::python::list params, ProjUnit units,
		       ProjDatum dat)
		{
			double p[16];
			for (int i=0; i<15; ++i)
				p[i] = extract<double>(params[i]);
			return;
		}
	ProjectionWrap(std::string params)
		{
			return;
		}
	ProjectionWrap(const Projection& s)
		{
			return;
		}
	void _init()
		{
			this->get_override("_init")();
			return;
		}
	void _forward(double lon, double lat)
		{
			this->get_override("_forward")();
			return;
		}	
	void _inverse(double lon, double lat)
		{
			this->get_override("_inverse")();
			return;
		}
};

BOOST_PYTHON_MODULE(pygctp)
{
	enum_<ProjCode>("ProjCode")
		.value("NONE", NONE)
		.value("GEO", GEO)
		.value("UTM", _UTM)
		.value("SPCS", SPCS)
		.value("ALBERS", ALBERS)
		.value("LAMCC", LAMCC)
		.value("MERCAT", MERCAT)
		.value("PS", PS)
		.value("POLYC", POLYC)
		.value("EQUIDC", EQUIDC)
		.value("TM", TM)
		.value("STEREO", STEREO)
		.value("LAMAZ", LAMAZ)
		.value("AZMEQD", AZMEQD)
		.value("GNOMON", GNOMON)
		.value("ORTHO", ORTHO)
		.value("GVNSP", GVNSP)
		.value("SNSOID", SNSOID)
		.value("EQRECT", EQRECT)
		.value("MILLER", MILLER)
		.value("VGRINT", VGRINT)
		.value("HOM", HOM)
		.value("ROBIN", ROBIN)
		.value("SOM", SOM)
		.value("ALASKA", ALASKA)
		.value("GOOD", GOOD)
		.value("MOLL", MOLL)
		.value("IMOLL", IMOLL)
		.value("HAMMER", HAMMER)
		.value("WAGIV", WAGIV)
		.value("WAGVII", WAGVII)
		.value("OBEQA", OBEQA)
		.value("USDEF", USDEF)
		;
	
	enum_<ProjUnit>("ProjUnit")
		.value("UNDEF", UNDEF)
		.value("RADIAN", RADIAN)
		.value("FEET", FEET)
		.value("METER", METER)
		.value("SECOND", SECOND)
		.value("DEGREE", DEGREE)
		.value("INT_FEET", INT_FEET)
		;

	enum_<ProjDatum>("ProjDatum")
		.value("NOT_SET", NOT_SET)
		.value("CLARKE_1866", CLARKE_1866)
		.value("BESSEL", BESSEL)
		.value("INTERNATIONAL_1967", INTERNATIONAL_1967)
		.value("INTERNATIONAL_1909", INTERNATIONAL_1909)
		.value("WGS_72", WGS_72)
		.value("EVEREST", EVEREST)
		.value("WGS_66", WGS_66)
		.value("GRS_1980", GRS_1980)
		.value("AIRY", AIRY)
		.value("MODIFIED_EVEREST", MODIFIED_EVEREST)
		.value("MODIFIED_AIRY", MODIFIED_AIRY)
		.value("WGS_84", WGS_84)
		.value("SOUTHEAST_ASIA", SOUTHEAST_ASIA)
		.value("AUSTRALIAN_NATIONAL", AUSTRALIAN_NATIONAL)
		.value("KRASSOVSKY", KRASSOVSKY)
		.value("HOUGH", HOUGH)
		.value("MERCURY_1960", MERCURY_1960)
		.value("MODIFIED_MERCURY_1968", MODIFIED_MERCURY_1968)
		.value("EARTH", EARTH)
		.value("BESSEL_1841_NAMIBIA", BESSEL_1841_NAMIBIA)
		.value("EVEREST_SABAH", EVEREST_SABAH)
		.value("EVEREST_INDIA_1956", EVEREST_INDIA_1956)
		.value("EVEREST_MALAYSIA_1969", EVEREST_MALAYSIA_1969)
		.value("EVEREST_MALAY_1948", EVEREST_MALAY_1948)
		.value("EVEREST_PAKISTAN", EVEREST_PAKISTAN)
		.value("HAYFORD", HAYFORD)
		.value("HELMERT_1906", HELMERT_1906)
		.value("INDONESIAN_1974", INDONESIAN_1974)
		.value("SOUTH_AMERICAN_1969", SOUTH_AMERICAN_1969)
		.value("WGS_60", WGS_60)
		;
	
	class_<ProjectionWrap>("Projection")
		.def(init<boost::python::list, ProjUnit, ProjDatum>())
		.def("forward", &ProjectionWrap::forward)
		.def("inverse", &Projection::inverse)
		.def("x", &Projection::x)
		.def("y", &Projection::y)
		.def("lat", &Projection::lat)
		.def("lon", &Projection::lon)

		.def("lat_lon", &Projection::latLon)
		.def("param", &Projection::param)
		.def("name", &Projection::name)
		.def("units", &Projection::units)
		.def("datum", &Projection::datum)
		.def("set_units", &Projection::setUnits)
		.def("set_datum", &Projection::setDatum)
		.def("set_fe", &Projection::setFE)
		.def("set_fn", &Projection::setFN)
		.def("set_radii", &Projection::setRadii)
//		.def{"set_major_radius", &Projection::setRMajor)
//		.def("set_minor_radius", &Projection::setRMinor)
//		.def("set_radius", &Projection::setRadius)
		;
	def("Projection", Transformer::convertProjection, "Returns a projection");

}
