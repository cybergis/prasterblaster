
#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "rasterinfo.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mercator.h"
#include "gctp_cpp/constants.h"
#include "gctp_cpp/utm.h"
#include "reprojector.hh"

#include "projectedraster.hh"

using namespace std;

ProjectedRaster::ProjectedRaster(string _filename)
{
	OGRSpatialReference sr;
	size_t found;

	projection = 0;
	filename = _filename;
	projection = 0;
	dataset = 0;
	data = 0;
	rows = cols = -1;

	GDALAllRegister();

	// Determine whether we have a imagine binary raster or a GTiff 
	found = filename.rfind("img");

	if (found == string::npos) { // Filename doesn't end in img, assuming GTiff
		ready = loadRaster(filename);
	} else { // Filename specifies imagine binary raster
		_filename.replace(found-1, 4, ".xml");
		bool status = configureFromXml(_filename);
		band_count = 1;
		ready = status;
		
	}
	

	return;
}

ProjectedRaster::ProjectedRaster(string _filename, 
				 int num_rows, int num_cols, 
				 GDALDataType pixel_type, double _pixel_size,
				 int _band_count,
				 Projection *proj,
				 double ulx, double uly)
{
	filename = _filename;
	projection = proj;
	band_count = _band_count;
	rows = num_rows;
	cols = num_cols;
	pixel_size = _pixel_size;
	type = pixel_type;
	ul_x = ulx;
	ul_y = uly;

	dataset = 0;

	ready = makeRaster();

	return;
}


ProjectedRaster::ProjectedRaster(ProjectedRaster *input,
				 string _filename,
				 string xmlDescription)
{
	projection = 0;
	filename = _filename;
	projection = 0;
	dataset = 0;
	data = 0;
	rows = cols = -1;

	GDALAllRegister();
	bool status = configureFromXml(xmlDescription);
	band_count = 1;

	if (status == true) {
		Area a = FindMinBox(input->ul_x, input->ul_y,
				    input->pixel_size,
				    rows, cols,
				    input->projection,
				    projection,
				    pixel_size);
	}
	

	ready = makeRaster();

	return;
}

ProjectedRaster::ProjectedRaster(string _filename,
				 ProjectedRaster *input,
				 Projection *output_proj,
				 GDALDataType pixel_type,
				 double _pixel_size)
{
	double lr_x, lr_y;
	const char *format = "GTiff";
	GDALDriver *driver;
	GDALRasterBand *band;
	Area area;

	projection = output_proj->copy();
	filename = _filename;
	dataset = 0;
	data = 0;
	band = 0;
	band_count = input->bandCount();
	rows = -1;
	cols = -1;
	pixel_size = _pixel_size;
	type = pixel_type;

	area = FindMinBox(input->ul_x, input->ul_y, input->pixel_size,
			  input->rows, input->cols, input->projection,
			  output_proj, pixel_size);

	ul_x = area.ul.x;
	ul_y = area.ul.y;
	lr_x = area.lr.x;
	lr_y = area.lr.y;

	rows = (int)ceil((ul_y-lr_y) / input->getPixelSize());
	cols = (int)ceil((lr_x-ul_x) / input->getPixelSize());

	GDALAllRegister();

	driver = GetGDALDriverManager()->GetDriverByName(format);
/*
	if( driver == NULL ) {
		ready = false;
		return;
	}
	
	if(filename != "") {
		// Set options
		options = CSLSetNameValue( options, "INTERLEAVE", "PIXEL" );
		options = CSLSetNameValue( options, "BIGTIFF", "NO" );
		options = CSLSetNameValue( options, "TILED", "NO" );
		options = CSLSetNameValue( options, "COMPRESS", "NONE" );
		options = CSLSetNameValue( options, "PHOTOMETRIC", "MINISWHITE");
		options = CSLSetNameValue( options, "PROFILE", "GDALGeoTiff");

		dataset = driver->Create(filename.c_str(), cols, rows,
					 band_count, pixel_type,
					 options);

	
		if (dataset == 0 || projection == 0) {
			ready = false;
			return;
		}
		
		geotransform[0] = ul_x;
		geotransform[3] = ul_y;
		geotransform[4] = geotransform[2] = 0.0;
		geotransform[1] = geotransform[5] = pixel_size;

		dataset->SetGeoTransform(geotransform);
	
		band = dataset->GetRasterBand(1);


		// Setup georeferencing
		OGRSpatialReference srs;
		long zone = -1;
		long projsys, datum;
		double *params;


		CPLErr err = dataset->SetProjection(projection->wkt().c_str());
		if (err = CE_None) {
			
		} else if (err = CE_Failure)  {
			fprintf(stderr, "Error setting projection."
			       "Was GDAL compiled with GTiff support?\n");
			ready = false;
		}
		dataset->FlushCache();
		GDALClose(dataset);
		CPLFree(wkt);

		dataset = (GDALDataset*)GDALOpen(filename.c_str(), GA_Update);


		ready = true;

		if (options != 0)
			CSLDestroy(options);
	}
*/
	
	ready = false;
	return;
}




ProjectedRaster::~ProjectedRaster()
{
	if (data != 0) {
		free(data);
		data = 0;
	}

	if (projection != 0) {
		free(projection);
		projection = 0;
	}

	if (dataset != 0) {
		GDALClose( (GDALDatasetH) dataset);
		dataset = 0;
	}


	return;
}

bool ProjectedRaster::isReady()
{
	if ((projection != 0) && (projection->error() == 0)
	    && ready) {
		return true;
	}  else {
		return false;
	}
}

Projection* ProjectedRaster::getProjection()
{
	return projection->copy();
}


int ProjectedRaster::getRowCount()
{
	return rows;
}

int ProjectedRaster::getColCount()
{
	return cols; 
}


GDALDataType ProjectedRaster::getPixelType()
{
	return type;
}

int ProjectedRaster::bitsPerPixel() 
{
	switch(getPixelType()) {
	case GDT_Byte:
		return 8;
	case GDT_Int16:
	case GDT_UInt16:
		return 16;
	case GDT_Int32:
	case GDT_UInt32:
	case GDT_Float32:
		return 32;
	case GDT_Float64:
		return 64;
	default:
		return -1;

	}


	return -1;
		
}

int ProjectedRaster::bandCount()
{
	if (dataset == 0) {
		return 1;
	} else {
		return dataset->GetRasterCount();
	}
	
	return -1;
	
}

double ProjectedRaster::getPixelSize()
{
	return pixel_size;
}

int ProjectedRaster::getZoneNumber()
{
	if (projection != 0) {
		if (projection->number() == _UTM) {
			return ((UTM*)projection)->zone();
		}
		
	}
	return -1;
}

ProjDatum ProjectedRaster::getDatum()
{
	if (projection != 0) {
		return projection->datum();
	}
	return (ProjDatum)-1;
}

double* ProjectedRaster::getGctpParams()
{
	return gctpParams;
}


bool ProjectedRaster::readRaster(int firstRow, int numRows, void *data)
{
	bool success = false;
	

	if ((firstRow + numRows) > getRowCount()) {
		printf("Brrrrrangn\n\n\n");
		return false;
	}

	if (isReady() && dataset == 0) { // img file
		ifstream ifs(filename.c_str(), ifstream::in);
		// TODO: Verify parameters!

		if (ifs.good()) {
			ifs.seekg(firstRow * cols * (GDALGetDataTypeSize(type)/8));
			ifs.read((char*)data, 
				 numRows * cols * (GDALGetDataTypeSize(type)/8));
		}
		if (ifs.good()) {
			success = true;
		} else if (ifs.fail()) {
			printf("File read failed!\n");
		} else if( ifs.bad()) {
			printf("Bad file read\n");
		} else if ( ifs.eof() ) {
			printf("Attempted read past end of file.\n");
		}
			
	} else if (isReady() && dataset != 0) { // GTiff
		if( dataset->RasterIO(GF_Read, 0, firstRow,
				      cols, numRows,
				      data, cols, numRows,
				      type, bandCount(),
				      NULL, 0, 0, 0) == CE_None) {
			success = true; 
		}
			
	}

		return success;
}

bool ProjectedRaster::writeRaster(int firstRow, int numRows, void *data)
{	
	bool success = false;

	if (firstRow < 0 || numRows +firstRow > rows || data == 0) {
		fprintf(stderr,"Write Boink #1\n");
		fflush(stderr);
			
		return false;
	}

	if (firstRow < 0 || firstRow >= dataset->GetRasterYSize()) {
		fprintf(stderr, "Write boink #2\n");
		fflush(stderr);
		return false;
	}

	if ((firstRow + numRows) > dataset->GetRasterYSize()) {
		fprintf(stderr, "Write boink #3\n");
		fflush(stderr);

		printf("Size: %d, Actual Size: %d\n", firstRow + numRows,
		       dataset->GetRasterYSize());
		return false;
	}
		
	
	if (isReady() && dataset == 0) { // img file
		ofstream ofs(filename.c_str(), ifstream::out);

		// TODO: Verify parameters!
		

		if (ofs.good()) {
			ofs.seekp(firstRow * cols);
			ofs.write((char*)data, 
				 numRows * cols * (GDALGetDataTypeSize(type)/8));
		} 
	} else if (isReady() && dataset != 0) { // GTiff
		if( dataset->RasterIO(GF_Write, 0, firstRow,
				      cols, numRows,
				      data, cols, numRows,
				      type, bandCount(),
				      NULL, 0, 0, 0) == CE_None) {
			success = true;
			dataset->FlushCache();
		}
			
	}

		return success;

	
	return false;
}

bool ProjectedRaster::configureFromXml(std::string xmlfilename)
{
	RasterInfo in_info(xmlfilename.c_str());

	if (!in_info.ready()) {
		return false;
	}


	rows = in_info.rows();
	cols = in_info.cols();
	
	ul_x = in_info.ul_X();
	ul_y = in_info.ul_Y();

	pixel_size = in_info.pixelSize();

	for(int i=0; i<15; ++i) {
		gctpParams[i] = in_info.gctpParam(i);
	}	

	if (in_info.dataType() == "Integer") {
		if (in_info.isSigned()) {
			switch (in_info.bitCount()) {
			case 8:
				type = GDT_Byte;
				break;
			case 16:
				type = GDT_Int16;
				break;
			case 32:
				type = GDT_Int32;
				break;
			default:
				return false;
			}
		} else { // Type is unsigned
			switch (in_info.bitCount()) {
			case 8:
				type = GDT_Byte;
				break;
			case 16:
				type = GDT_UInt16;
				break;
			case 32:
				type = GDT_UInt32;
				break;
			default: 
				return false;
			}
		}
	} else { // Type is float
		if (in_info.bitCount() == 32)
			type = GDT_Float32;
		else
			type = GDT_Float64;
	}

	// Setup projection
	projection = 
	  Transformer::convertProjection((ProjCode)in_info.projectionNumber());

	if (projection == 0) {

		return false;
	}
	projection->setUnits((ProjUnit)in_info.unitNumber());
	projection->setDatum((ProjDatum)in_info.datumNumber());
	projection->setParams(in_info.allGctpParams());

	return true;
}

bool ProjectedRaster::loadImgRaster(string rasterFilename, string xmlFilename)
{
	rasterFilename = "";
	xmlFilename  = "";
	
	return false;
}

bool ProjectedRaster::loadRaster(string filename)
{
	projection = 0;
	projection = 0;
	dataset = 0;
	data = 0;
	rows = cols = -1;
	char *ref, **ugh;
	long projsys, zone, datum;
	double params[18];
	double *p;
	OGRSpatialReference sr;
	

	GDALAllRegister();	
	
	dataset = (GDALDataset*)GDALOpen( filename.c_str(), GA_Update );
	
	if (dataset == 0) {
		return false;
	} 
	
	cols = dataset->GetRasterXSize();
	rows = dataset->GetRasterYSize();
	
	
	ref = strdup(dataset->GetProjectionRef());
	
	if (ref == 0) {
		fprintf(stderr, "Error reading projection ref\n");
		return false;
	}
	
	double geo[6];
	dataset->GetGeoTransform(geo);
	pixel_size = geo[1];
	ul_x = geo[0];
	ul_y = geo[3];
	type = (dataset->GetRasterBand(1))->GetRasterDataType();
	band_count = dataset->GetRasterCount();
	cols = dataset->GetRasterXSize();
	rows = dataset->GetRasterYSize();

	// Setup projection
	ugh = &ref;
	sr.importFromWkt(ugh);
	
	p = &(params[0]);
	sr.exportToUSGS(&projsys, &zone, &p, &datum);
	projection = Transformer::convertProjection((ProjCode)projsys);
	if (projection == 0) {
		fprintf(stderr, "Error building projection, num %ld...\n", projsys);
		return false;
	}
	if ((ProjCode)projsys == _UTM)
		((UTM*)projection)->setZone(zone);
	projection->setParams(params);
	projection->setDatum((ProjDatum)datum);

	projection->setUnits(METER);
	
	return true;
}

bool ProjectedRaster::makeRaster()
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
	options = CSLSetNameValue( options, "BIGTIFF", "NO" );
	options = CSLSetNameValue( options, "TILED", "NO" );
	options = CSLSetNameValue( options, "COMPRESS", "NONE" );
	options = CSLSetNameValue( options, "PHOTOMETRIC", "MINISBLACK");

	dataset = driver->Create(filename.c_str(), cols, rows, 
				 band_count, type,
				 options);

	
	if (dataset == 0 || projection == 0) {
		return false;
	}

	// Setup georeferencing
	OGRSpatialReference srs;

	geotransform[0] = ul_x;
	geotransform[3] = ul_y;
	geotransform[4] = geotransform[2] = 0.0;
	geotransform[1] = geotransform[5] = pixel_size;
	dataset->SetGeoTransform(geotransform);
	dataset->SetProjection(projection->wkt().c_str());


	if (options != 0)
		CSLDestroy(options);

	return true;

}
