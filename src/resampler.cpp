

#include "resampler.hh"

namespace resampler
{

// This function is to coax the compiler to actually generate some
// code. Don't call it.
  void never_call_me()
  {
	  nearest_neighbor<unsigned char>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<signed char>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<int>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<unsigned int>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<long>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<unsigned long>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<float>(0,0,0,0,0,0,0,0);
	  nearest_neighbor<double>(0,0,0,0,0,0,0,0);
  }
  
}
