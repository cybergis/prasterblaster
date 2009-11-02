

#ifndef RESAMPLER_HH
#define RESAMPLER_HH

template <typename T>
void nearest_neighbor(void *inraster, double in_x,
		      double in_y, unsigned long in_cols,
		      void *outraster,
		      unsigned long out_x, unsigned long out_y,
		      unsigned long out_cols);

		      


#endif // RESAMPLER_HH
