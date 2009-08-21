

#include <boost/python.hpp>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h" // To include all subclasses of Projection

using namespace boost::python;

struct ProjectionWrap : Projection, wrapper<Projection>
{
	ProjectionWrap() : Projection() {};
	ProjectionWrap(double args[], ProjUnit units, ProjDatum dat) : Projection(args, units, dat) {};
	void _init()
		{
			this->get_override("_init")();
		}
	void _forward(double lon, double lat)
		{
			this->get_override("_forward")(lon, lat);
		}
	void _inverse(double lon, double lat)
		{
			this->get_override("_inverse")(lon, lat);
		}

};

ProjCode identity_(ProjCode x) { return x; }


BOOST_PYTHON_MODULE(Projection)
{
	enum_<ProjCode>("ProjCode")
		.value("NONE", NONE)
		.value("GEO", GEO)
		.value("_UTM", _UTM)
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
		;
	def("identity", identity_);

	class_<ProjectionWrap, boost::noncopyable>("Projection")
//		.def(init<double[15], int, int>)
//		.def("_init", pure_virtual(&Projection::_init))
//		.def("_forward", pure_virtual(&Projection::_forward))
//		.def("_inverse", pure_virtual(&Projection::_inverse))
		.def("forward", &Projection::forward)
		.def("inverse", &Projection::inverse)
		;
	class_<Equirectangular, bases<Projection> >("Equirectangular")
		;
}
