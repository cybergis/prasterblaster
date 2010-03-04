

#ifndef RESAMPLER_HH
#define RESAMPLER_HH

#include <gdal/gdal.h>

namespace resampler 
{
	enum RESAMPLERS {
		NO_RESAMPLER = 0,
		NEAREST_NEIGHBOR,
		BILINEAR
	};

	typedef void(*resampler_func)(void *inraster, double in_x,
				      double in_y, unsigned long in_cols,
				      void *outraster,
				      unsigned long out_x, unsigned long out_y,
				      unsigned long out_cols);

	resampler_func SelectResampler333(RESAMPLERS method, GDALDataType type);

		
	template <typename T>
	void nearest_neighbor(void *inraster, double in_x,
			      double in_y, unsigned long in_cols,
			      void *outraster,
			      unsigned long out_x, unsigned long out_y,
			      unsigned long out_cols);
	
	
	template <typename T>
	void bilinear(void *inraster, double in_x,
		      double in_y, unsigned long in_cols,
		      void *outraster,
		      unsigned long out_x, unsigned long out_y,
		      unsigned long out_cols);
	
} // namespace resampler

#endif // RESAMPLER_HH
