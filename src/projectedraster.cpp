/*!
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE
 * This work was produced as a part of the official duties of a
 * federal employee and is in the public domain.

 * @section DESCRIPTION
 *
 * The ProjectedRaster class represents a raster with a location and a projection.
 *
 */

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
using std::shared_ptr;


ProjectedRaster::ProjectedRaster(string _filename)
{
	OGRSpatialReference sr;
	size_t found;

	filename = _filename;
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

bool ProjectedRaster::CreateRaster(string _filename, 
				   int num_rows, int num_cols, 
				   GDALDataType pixel_type, double _pixel_size,
				   int _band_count,
				   shared_ptr<Projection> proj,
				   double ulx, double uly)
{
	bool status = makeRaster(_filename,
				 num_cols,
				 num_rows,
				 _band_count,
				 ulx,
				 uly,
				 pixel_type,
				 proj,
				 _pixel_size);
			   
	return status;
}


bool ProjectedRaster::CreateRaster(shared_ptr<ProjectedRaster> input,
				   string _filename,
				   string xmlDescription)
{

	RasterInfo in_info(xmlDescription.c_str());
	double gctpParams[16];
	GDALDataType type;
	Projection *projection = 0;
	int _band_count = 1;

	input->getColCount();
	if (!in_info.ready()) {
		fprintf(stderr, "Could not open xml file: %s\n",
			xmlDescription.c_str());
		return false;
	}


	int num_rows = in_info.rows();
	int num_cols = in_info.cols();
	double ulx = in_info.ul_X();
	double uly = in_info.ul_Y();

	double _pixel_size = in_info.pixelSize();

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
		fprintf(stderr, "Error creating projection for raster: %s\n",
			_filename.c_str());
		return false;
	}
	projection->setUnits((ProjUnit)in_info.unitNumber());
	projection->setDatum((ProjDatum)in_info.datumNumber());
	projection->setParams(in_info.allGctpParams());

	bool status = makeRaster(_filename,
				 num_cols,
				 num_rows,
				 _band_count,
				 ulx,
				 uly,
				 type,
				 shared_ptr<Projection>(projection->copy()),
				 _pixel_size);

			  

	return status;
}


bool ProjectedRaster::CreateRaster(string _filename,
				   shared_ptr<ProjectedRaster> input,
				   shared_ptr<Projection> output_proj,
				   GDALDataType pixel_type,
				   double _pixel_size)
{
	bool status = false;
	double ulx, uly;
	int num_cols, num_rows;

	Area out_area =  FindOutputArea(input,
					  output_proj,
					  _pixel_size);
	ulx = out_area.ul.x;
	uly = out_area.ul.y;

	num_rows = (int)((out_area.ul.y - out_area.lr.y) / _pixel_size);
	num_cols = (int)((out_area.lr.x - out_area.ul.x) / _pixel_size);
	printf("size is %d cols %d rows\n", num_cols, num_rows);

	status = makeRaster(_filename,
			    num_cols, 
			    num_rows,
			    input->band_count,
			    ulx,
			    uly,
			    pixel_type,
			    output_proj,
			    _pixel_size);



	return status;
}




ProjectedRaster::~ProjectedRaster()
{
	if (data != 0) {
		free(data);
		data = 0;
	}

	if (dataset != 0) {
		GDALClose( (GDALDatasetH) dataset);
		dataset = 0;
	}


	return;
}

bool ProjectedRaster::isReady()
{
  if (projection.get() != 0 && 
      (projection->error() == 0) && ready) {
		return true;
	}  else {
		return false;
	}
}

shared_ptr<Projection> ProjectedRaster::getProjection()
{
	return projection;
}


int ProjectedRaster::getRowCount()
{
	return rows;
}

int ProjectedRaster::getColCount()
{
	return cols; 
}

void ProjectedRaster::clampGeoCoordinate(Coordinate *c)
{
	if (c == 0) {
		return; 
	}

	if (c->x < -180.0) {
		c->x = -179.999999;
	}
	
	if (c->x > 180.0) {
		c->x = 179.9999999;
	}

	if (c->y < -90.0) {
		c->y = -89.9999999;
	}


	if (c->y > 90.0) {
		c->y = 89.99999999;
	}

	return;
}

Coordinate ProjectedRaster::getGeographicalCoordinate(int rasterX, int rasterY)
{
	Coordinate c = getProjectedCoordinate(rasterX, rasterY);

	projection->inverse(c.x, c.y, &c.x, &c.y);

	return c;
	
}

Coordinate ProjectedRaster::getProjectedCoordinate(int rasterX, int rasterY)
{
	Coordinate c;

	c.x = ul_x + rasterX * pixel_size;
	c.y = ul_y - rasterY * pixel_size;

	return c;

}

int ProjectedRaster::getPixelIndex(double longitude, double latitude)
{

	double x, y;
	projection->forward(longitude, latitude, &x, &y);
	x = ul_x - x;
	y = ul_x - y;
	latitude += ul_y;
	longitude /= pixel_size;
	latitude /= pixel_size;

	return 0;

}

int ProjectedRaster::getPixelIndex(Coordinate geographicalCoordinate)
{

	return ProjectedRaster::getPixelIndex(geographicalCoordinate.x, 
					      geographicalCoordinate.y);
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
/*
	shared_ptr<UTM> utm_projection(projection);
	
	if (projection != 0) {
		if (projection->number() == _UTM) {
			return (utm_projection->zone());
		}
		
	}
*/
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
			return false;
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
	  fprintf(stderr,"Write Boink #1 %d %d\n", numRows +firstRow, rows);
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
		shared_ptr<Projection>(Transformer::convertProjection((ProjCode)in_info.projectionNumber()));

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
	dataset = 0;
	data = 0;
	rows = cols = -1;
	char *ref, **ugh;
	long projsys, zone, datum;
	double params[18] = {0.0};
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
	projection = shared_ptr<Projection>(Transformer::convertProjection((ProjCode)projsys));
	if (projection == 0) {
		fprintf(stderr, "Error building projection, num %ld...\n", projsys);
		return false;
	}
/*	shared_ptr<UTM> utm_projection(projection);
	if ((ProjCode)projsys == _UTM)
	(utm_projection->setZone(zone); */
	projection->setParams(params);
	projection->setDatum((ProjDatum)datum);

	projection->setUnits(METER);
	ready = true;
	
	return true;
}

bool ProjectedRaster::makeRaster(string _filename,
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
	options = CSLSetNameValue( options, "BIGTIFF", "YES" );
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

	GDALClose(_dataset);

	if (options != 0)
		CSLDestroy(options);

	return true;

}
