/* 
 * Programmer: David Mattli <dmattli@usgs.gov>
 */


#include <string>

using std::string;

int driver(string input_raster, 
	   string output_filename, 
	   string temporary_path,
	   string output_srs, 
	   string resampler,
	   string fillvalue,
	   int partition_count);
