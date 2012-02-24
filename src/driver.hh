/* 
 * Programmer: David Mattli <dmattli@usgs.gov>
 */


#include <string>

using std::string;

int driver(string input_raster, 
	   string output_filename, 
	   string output_srs, 
	   string fillvalue,
	   int partition_count);
