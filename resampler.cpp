

#include "resampler.hh"

template <typename type>
void nearest_neighbor(void *inraster, double in_x,
		      double in_y, unsigned long in_cols,
		      void *outraster,
		      unsigned long out_x, unsigned long out_y,
		      unsigned long out_cols)
{
  
  ((type*)outraster)[out_x + (out_y * out_cols)] = 
    ((type*)inraster)[(unsigned long)in_x + ((unsigned long)in_y * in_cols)];
  return;
}


// This function is to coax the compiler to actually generate some
// code. Don't call it.
void never_call_me()
{
  nearest_neighbor<unsigned char>(0,0,0,0,0,0,0,0);
  nearest_neighbor<signed char>(0,0,0,0,0,0,0,0);
  nearest_neighbor<int>(0,0,0,0,0,0,0,0);
  nearest_neighbor<unsigned int>(0,0,0,0,0,0,0,0);
  nearest_neighbor<float>(0,0,0,0,0,0,0,0);
  nearest_neighbor<double>(0,0,0,0,0,0,0,0);
}
