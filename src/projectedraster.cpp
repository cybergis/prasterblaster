
#include <QString>

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
	char *ref = 0;
	char **ugh = 0;
	double params[16], *first;
	long int zone, projsys, datum;
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

		dataset = (GDALDataset*)GDALOpen( _filename.c_str(), GA_ReadOnly );

		if (dataset == 0) {
			ready = false;
			return;
		} 
	
		ref = strdup(dataset->GetProjectionRef());

		if (ref == 0) {
			ready = false;
			return;
		}

		double geo[6];
		geo[1] = -1;
		dataset->GetGeoTransform(geo);
		pixel_size = geo[1];
		
		// Setup projection
		ugh = &ref;
		sr.importFromWkt(ugh);
		first = &(params[0]);
		sr.exportToUSGS(&projsys, &zone, &first, &datum);
		projection = Transformer::convertProjection((ProjCode)projsys);
		if ((ProjCode)projsys == _UTM)
			((UTM*)projection)->setZone(zone);
		projection->setParams(params);
		projection->setDatum((ProjDatum)datum);
		projection->setUnits(METER);

		free(ref);
		ready = true;
		
	} else { // Filename specifies imagine binary raster
		bool status = loadImgRaster(filename);
		ready = status;
	}
	
	return;
}

ProjectedRaster::ProjectedRaster(string _filename, 
				 int num_rows, int num_cols, 
				 GDALDataType pixel_type, double pixel_size,
				 int band_count,
				 Projection *proj,
				 double ulx, double uly)
{
	const char *format = "GTiff";
	GDALDriver *driver;
	GDALRasterBand *band;
	char **options = 0; 
	GByte raster[8*8];
	char *wkt = 0;
	double geotransform[6] = { 444720, 30, 0, 3751320, 0, -30 };

	projection = 0;
	filename = _filename;
	projection = 0;
	dataset = 0;
	data = 0;
	band = 0;
	rows = cols = -1;

	GDALAllRegister();

	projection = proj;

	driver = GetGDALDriverManager()->GetDriverByName(format);

	if( driver == NULL ) {
		ready = false;
		return;
	}
	  
	// Set options
//	options = CSLSetNameValue( options, "INTERLEAVE", "PIXEL" );
//	options = CSLSetNameValue( options, "BIGTIFF", "YES" );
//	options = CSLSetNameValue( options, "TILED", "YES" );
//	options = CSLSetNameValue( options, "COMPRESS", "NONE" );

	printf("Creating dataset! %d cols, %d rows, %d bands\n",
	       num_cols, num_rows, band_count);
	dataset = driver->Create(filename.c_str(), num_cols, num_rows, 
				 band_count, pixel_type,
				 0);
/*
	if (dataset == 0 || proj == 0) {
		ready = false;
		return;
	}

	dataset->SetGeoTransform(geotransform);
	
	band = dataset->GetRasterBand(1);

	band->RasterIO(GF_Write, 0, 0, 8, 8, raster, 8, 8, GDT_Byte, 0, 0);

	// Setup georeferencing
	OGRSpatialReference srs;
	long zone = -1;
	
	if (projection->number() == _UTM) {
		zone = (int)((UTM*)projection)->zone();
	}
	srs.importFromUSGS((long)projection->number(), 
			    zone,
			    projection->params(),
			    (long)projection->datum(),
			    TRUE);
	srs.exportToWkt(&wkt);
	printf("WKT YO: %s\n", wkt);
	dataset->SetProjection(wkt);
	CPLFree(wkt);
		
	

	ready = true;

	if (options != 0)
		CSLDestroy(options);
*/
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
	if ((projection != 0) && (projection->error() == 0)) {
		return true;
	}  else {
		return false;
	}
}

Projection* ProjectedRaster::getProjection()
{
	return projection->copy();
}

void* ProjectedRaster::getData()
{
	return data;
}

bool ProjectedRaster::write(string filename)
{
	GDALDataset *ds;
	GDALDriver *driver;
	OGRSpatialReference sr;
	Projection *proj;
	char **options = 0;
	char *wkt;
	double transforms[6];

	GDALAllRegister();

	driver = GetGDALDriverManager()->GetDriverByName("GTiff");
/*  
	// Initialize spatial reference
	// Projection, gctp params, units, and datum
	proj = getProjection();
	if (getProjection()->number() == _UTM) {
		sr.importFromUSGS((long)proj->number(), ((UTM*)proj)->zone(),
				  proj->params(), (long)proj->datum());
	} else {
		sr.importFromUSGS((long)proj->number(),  0,
				  proj->params(), (long)proj->datum());
	}
	sr.SetLinearUnits(SRS_UL_METER, 1);
	sr.exportToWkt( &wkt );

	// Create Dataset
	ds = driver->Create(filename.c_str(), getColCount(),
			    getRowCount(), getBitCount()/8, 
			    GDT_Byte, options);
  


	if (ds != 0) {
		// TODO: Support more than one band
		printf("Writing raster: %d cols %d rows\n", getColCount(), getRowCount());
		if (ds->RasterIO(GF_Write, 0, 0, getColCount(), getRowCount(),
				 getData(), getColCount(), getRowCount(),
				 GDT_Byte, 1, 0, 0, 0, 0) == CE_Failure) {
			printf("RasterIO failed...\n");
			return false;
		}

		ds->SetProjection(wkt);
		GDALClose(ds);
		return true;
	}
*/	
	return false;

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
	GDALRasterBand *band = 0;
	GDALDataType type = GDT_Unknown;
	
	if ( dataset == 0) {
		return type;
	}

	band = dataset->GetRasterBand(1);
	if (band == 0) {
		return type;
	} 
	
	return band->GetRasterDataType();
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


bool ProjectedRaster::loadImgRaster(std::string filename)
{
	string::size_type idx = filename.find('.');
	string imgname = filename.substr(0, idx) + ".img";
	string xmlname = filename.substr(0, idx) + ".xml";
	RasterInfo in_info(xmlname.c_str());
	int rastersize = in_info.rows() * in_info.cols();
	
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
	projection = Transformer::convertProjection((ProjCode)in_info.projectionNumber());
	if (projection == 0) {
		return false;
	}
	projection->setUnits((ProjUnit)in_info.unitNumber());
	projection->setDatum((ProjDatum)in_info.datumNumber());
	projection->setParams(in_info.allGctpParams());


	return true;
}

bool ProjectedRaster::readImgRaster(std::string filename)
{
	string::size_type idx = filename.find('.');
	string imgname = filename.substr(0, idx) + ".img";
	string xmlname = filename.substr(0, idx) + ".xml";
	RasterInfo in_info(xmlname.c_str());
	int fd = -1;
	int readsize = 0;
	int rastersize = in_info.rows() * in_info.cols();
	GDALDriver *driver = 0;
	char **options = 0;
	std::stringstream out;

	driver = GetGDALDriverManager()->GetDriverByName("MEM");

	if (driver == 0) {
		printf("Memory driver not found!\n");
		return false;
	}
	
	// Allocate memory for raster
	

	// Read in image data
	errno = 0;
	fd = open(imgname.c_str(), O_RDONLY);
	readsize = read(fd, data, 
			rastersize*(in_info.bitCount()/8));
	if (readsize != (rastersize * (in_info.bitCount()/8))) {
		printf("Read error: %d!\n", readsize);
		printf("Supposed to be: %d\n", rastersize*(in_info.bitCount()/8));
		printf("ERROR: %s\n", strerror(errno));
				return false;
	}
	close(fd);
	
	printf("Index: %d\n", idx);
	printf("image name: %s, xml name %s\n", imgname.c_str(), xmlname.c_str());

	out << "MEM:::DATAPOINTER=" << data
	    << ",PIXELS=" << in_info.cols()
	    << ",LINES=" << in_info.rows()
	    << ",BANDS=1" << "DATATYPE=Byte";
	dataset = driver->Create(out.str().c_str(), in_info.cols(), in_info.rows(), 
				 1, GDT_Byte, 0);

/*
	if (subIndex < subCount - 1)
		rastersize = subrastersize * in_info.cols();
	else
		rastersize = (subrastersize + suboverflow) * in_info.cols();
	
	data = calloc(rastersize, in_info.bitCount()/8);
	if (data == 0) {
		printf("Data allocation failed\n");
		return false;
	}

	// Read in image data
	errno = 0;
	fd = open(imgname.c_str(), O_RDONLY);
	lseek(fd, rastersize * subIndex, SEEK_SET);
	readsize = read(fd, data, 
			rastersize*(in_info.bitCount()/8));
	if (readsize != (rastersize * (in_info.bitCount()/8))) {
		printf("Read error: %d!\n", readsize);
		printf("Supposed to be: %d\n", rastersize*(in_info.bitCount()/8));
		printf("ERROR: %s\n", strerror(errno));
				return false;
	}
	close(fd);

	// Setup projection
	int num_rows = in_info.rows();
	

	setProjection((ProjCode)in_info.projectionNumber());
	setUL(in_info.ul_X(), in_info.ul_Y());
	setUL(in_info.ul_X(), 
	      in_info.ul_Y() - (subrastersize * in_info.pixelSize()));
	setDatum((ProjDatum)in_info.datumNumber());
	setRowCount(num_rows);
	if (subIndex < subCount -1) {
		setRowCount(subrastersize);
	} else {
		setRowCount(subrastersize + suboverflow);
	}
	setColCount(in_info.cols());
	setPixelSize(in_info.pixelSize());
	if (in_info.isSigned())
		setSigned();
	else
		setUnsigned();
	setUnit((ProjUnit)in_info.unitNumber());
	setGctpParams(in_info.allGctpParams());
	setBitCount(in_info.bitCount());
  
	return true;
*/
}

