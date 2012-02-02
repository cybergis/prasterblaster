

#include <cstdio>
#include <memory>
#include <string>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include <fcntl.h>

#include "projection.h"
#include "rasterinfo.h"
#include "transformer.h"

using std::string;
using std::shared_ptr;

bool makeRaster(string _filename,
		int _cols,
		int _rows,
		int _band_count,
		double _ul_x,
		double _ul_y,
		GDALDataType _type,
		shared_ptr<Projection> _projection,
		double _pixel_size)
				 
{
	const char *format = "GTiff";
	GDALDriver *driver;
	char **options = 0; 
	double geotransform[6] = { 444720, 30, 0, 3751320, 0, -30 };

	GDALAllRegister();

	driver = GetGDALDriverManager()->GetDriverByName(format);

	if( driver == NULL ) {
		return false;
	}
	  
	// Set options
	options = CSLSetNameValue( options, "INTERLEAVE", "PIXEL" );
//	options = CSLSetNameValue( options, "BIGTIFF", "YES" );
	options = CSLSetNameValue( options, "TILED", "NO" );
	options = CSLSetNameValue( options, "COMPRESS", "NONE" );
	options = CSLSetNameValue( options, "PHOTOMETRIC", "MINISBLACK");

	GDALDataset *_dataset = driver->Create(_filename.c_str(), _cols, _rows, 
				 _band_count, _type,
				 options);

	if (_dataset == 0 || _projection == 0) {
		return false;
	}

	// Setup georeferencing
	OGRSpatialReference srs;

	geotransform[0] = _ul_x;
	geotransform[3] = _ul_y;
	geotransform[4] = geotransform[2] = 0.0;
	geotransform[1] = geotransform[5] = _pixel_size;
	_dataset->SetGeoTransform(geotransform);

	_dataset->SetProjection(_projection->wkt().c_str());
	if (_projection->wkt() == "") {
		fprintf(stderr, "\nERROR: Unsupported projection: %s\n", _projection->name().c_str());

	}

	GDALClose(_dataset);

	if (options != 0)
		CSLDestroy(options);

	return true;

}

bool img2tif(char* input_filename, char* output_filename)
{
	size_t buffer_size = 0;
	char *buffer = NULL;
	ssize_t bytes_read = 0;
	ssize_t read_so_far = 0;
	int fd = open(input_filename, O_RDONLY);
	GDALDataset *ds = NULL;
	string xml_filename = input_filename;
	shared_ptr<Projection> projection; 

	if (fd < 0) {
		fprintf(stderr, "Error opening input file!\n");
		return false;
	}
	
	xml_filename.replace(xml_filename.rfind("img"), 4, "xml");

	RasterInfo in_info(xml_filename.c_str());

	if (!in_info.ready()) {
		fprintf(stderr, "Error reading xml file!\n");
		return false;
	}

	// Allocate buffer
	buffer_size = in_info.cols();
	buffer = (char*)malloc(buffer_size);

	projection = 
		shared_ptr<Projection>(Transformer::convertProjection((ProjCode)in_info.projectionNumber()));
	
	if (projection == 0) {
		return false;
	}
	projection->setUnits((ProjUnit)in_info.unitNumber());
	projection->setDatum((ProjDatum)in_info.datumNumber());
	projection->setParams(in_info.allGctpParams());

	if (makeRaster(output_filename, 
		       in_info.cols(),
		       in_info.rows(),
		       1,
		       in_info.ul_X(),
		       in_info.ul_Y(),
		       GDT_Byte,
		       projection,
		       in_info.pixelSize()) == false) {
		fprintf(stderr, "Error creating output raster!\n");
		return false;
	}
	
	ds = (GDALDataset*)GDALOpen(output_filename, GA_Update);

	if (ds == NULL) {
		return false;
	}

	bool done = false;
	while (!done) {
		bytes_read = read(fd, buffer, buffer_size);
		read_so_far += bytes_read;

		if (bytes_read > 0 ) {
			ds->RasterIO(GF_Write,
				     0,
				     read_so_far / in_info.cols(),
				     bytes_read,//				     read_so_far % in_info.cols(),
				     1,
				     buffer,
				     in_info.cols(),
				     1,
				     GDT_Byte,
				     1,
				     NULL, 0, 0, 0);
		} else {
			done = true;
			ds->FlushCache();
		}


	}
	

	// Clean up 
	close(fd);
	GDALClose(ds);


	return true;
}


int main(int argc, char *argv[]) 
{
	string usage = "USAGE: img2tif <input filename> <output filename>\n";

	if (argc < 3) {
		printf("%s\n", usage.c_str());
		return 0;
	}

	bool result = img2tif(argv[1], argv[2]);

	if (result == false) {
		fprintf(stderr, "Error occured!\n");
		return 1;
	}


	return 0;
}
