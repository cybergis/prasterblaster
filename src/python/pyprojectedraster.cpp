

#include <boost/python.hpp>

#include "../projectedraster.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(pRasterBlaster)
{
	class_<ProjectedRaster>("ProjectedRaster", init<std::string>())
		.def("isReady", &ProjectedRaster::isReady)
		.def("bandCount", &ProjectedRaster::bandCount)
		;
}
