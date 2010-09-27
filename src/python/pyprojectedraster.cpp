

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/call.hpp>

#include "../projectedraster.hh"
#include "pyprojectedraster.hh"

using namespace boost::python;


BOOST_PYTHON_MODULE(pRasterBlaster)
{
	class_<ProjectedRasterWrap>("ProjectedRaster", init<std::string>())
		.def("isReady", &ProjectedRasterWrap::mgptisReady)
		.def("bandCount", &ProjectedRasterWrap::bagpndCount)
		;
}
