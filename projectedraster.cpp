
#include <QString>

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "rasterinfo.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mercator.h"
#include "reprojector.hh"

#include "projectedraster.hh"

using namespace std;

ProjectedRaster::ProjectedRaster(string filename) 
{
	projection = 0;
	data = 0;
	if (readImgRaster(filename) == true) {
		ready = true;
	} else {
		ready = false;
	}
}

ProjectedRaster::ProjectedRaster(long num_rows, long num_cols, 
				 long pixel_bits, Projection *proj,
				 double ulx, double uly)
{
	printf("Allocatin' %ld\n", pixel_bits);
	data = calloc(num_cols*num_rows, (pixel_bits/8));
	ul_x = ulx;
	ul_y = uly;
	rows = num_rows;
	cols = num_cols;
	bitCount = pixel_bits;
	issigned = false;
	integer = true;
	projection = proj;

	for(int i = 0; i < 16; ++i) {
		gctpParams[i] = 0;
	}
	ready = true;

	if (data == 0)
		printf("Oh no!\n");

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

	return;
}

bool ProjectedRaster::isReady()
{
	if ((projection != 0) && (projection->error() == 0)
	    && (data != 0)) {
		return true;
	}  else {
		return false;
	}
}

Projection* ProjectedRaster::getProjection()
{
	return projection;
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
		}

		ds->SetProjection(wkt);
		GDALClose(ds);
		return true;
	}
	
	return false;

}

void ProjectedRaster::setUL(double _ul_x, double _ul_y)
{
	ul_x = _ul_x;
	ul_y = _ul_y;
	return;
}

void ProjectedRaster::setRowCount(int _rows)
{
	rows = _rows;
	return;
}

void ProjectedRaster::setColCount(int _cols)
{
	cols = _cols;
	return;
}

int ProjectedRaster::getRowCount()
{
	return rows;
}

int ProjectedRaster::getColCount()
{
	return cols;
}

void ProjectedRaster::setDataType(const std::string &datatype)
{
	QString str(datatype.c_str());
	type = datatype;
	if (str.contains("float") || str.contains("Float") 
	    || str.contains("IEEE"))
		setFloat();
	if (str.contains("8"))
		setBitCount(8);
	if (str.contains("32"))
		setBitCount(32);
	if (str.contains("unsigned") || str.contains("Unsigned"))
		setUnsigned();

	return;
}

string ProjectedRaster::getDataType()
{
	string t = "";

	if (bitCount == 8)
		t.append(" 8");
	else 
		t.append(" 32");

	if (isFloat())
		t.append("IEEE Float");
	else 
		t.append("Integer");

	if (!isSigned())
		t.append("Unsigned ");

	return t;
}

void ProjectedRaster::setPixelSize(double _pixsize)
{
	pixsize = _pixsize;
	return;
}

double ProjectedRaster::getPixelSize()
{
	return pixsize;
}

int ProjectedRaster::getBitCount() 
{
	return bitCount;
}


void ProjectedRaster::setBitCount(int _bits)
{
	bitCount = _bits;
	return;
}

void ProjectedRaster::setSigned()
{
	issigned = true;
	return;
}

void ProjectedRaster::setUnsigned()
{
	issigned = false;
	return;
}

bool ProjectedRaster::isSigned() 
{
	return issigned;
}

void ProjectedRaster::setInteger()
{
	integer = true;
	return;
}

bool ProjectedRaster::isInteger()
{
	if (integer == true)
		return true;
	else 
		return false;
}

void ProjectedRaster::setFloat()
{
	integer = false;
	return;
}

bool ProjectedRaster::isFloat()
{
	return !integer;
}

void ProjectedRaster::setProjection(ProjCode p)
{
	if (projection != 0) {
		delete projection;
		projection = 0;
	}
	projection = t.convertProjection(p);
	
	return;
}

void ProjectedRaster::setZoneNumber(int zone)
{
	if (projection != 0) {
		//    projection->setZone(zone);
	}
	return;
}

void ProjectedRaster::setUnit(ProjUnit unit) 
{
	if (projection != 0) {
		projection->setUnits(unit);
	}
	return;
}

void ProjectedRaster::setDatum(ProjDatum pd)
{
	if (projection != 0) {
		projection->setDatum(pd);
	}
	return;
}

void ProjectedRaster::setGctpParams(double params[])
{
	if (projection != 0) {
		projection->setParams(params);
	}
	return;
}

int ProjectedRaster::getZoneNumber()
{
	if (projection != 0) {
		//    Zone();
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

bool ProjectedRaster::readImgRaster(std::string filename)
{
	string::size_type idx = filename.find('.');
	string imgname = filename.substr(0, idx) + ".img";
	string xmlname = filename.substr(0, idx) + ".xml";
	RasterInfo in_info(xmlname.c_str());
	int fd = -1;
	int readsize = 0;
	int rastersize = in_info.rows() * in_info.cols();

	data = calloc(rastersize, in_info.bitCount()/8);
	if (data == 0) {
		printf("Data allocation failed\n");
		return false;
	}

	// Read in image data
	errno = 0;
	fd = open(imgname.c_str(), O_RDONLY);
	readsize = read(fd, data, 
			rastersize*(in_info.bitCount()/8));
	if (readsize != (rastersize * (in_info.bitCount()/8))) {
		printf("Read error: %d!\n", readsize);
		printf("Supposed to be: %d\n", rastersize*(in_info.bitCount()/8));
		printf("ERROR: %s\n", strerror(errno));
		return 0;
	}
	close(fd);

	// Setup projection
	int num_rows = in_info.rows();
	setProjection((ProjCode)in_info.projectionNumber());
	setUL(in_info.ul_X(), in_info.ul_Y());
	setDatum((ProjDatum)in_info.datumNumber());
	setRowCount(num_rows);
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
}
