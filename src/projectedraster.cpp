
#include <QString>

#include <string>
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

	projection = 0;
	filename = _filename;
	projection = 0;
	
	GDALAllRegister();

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

	// Else
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


        // clean up
	//free(ref);
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
	char **papszMetadata;
	char **options = 0;

	GDALAllRegister();

	projection = proj;
	filename = _filename;

	driver = GetGDALDriverManager()->GetDriverByName("GTiff");

	if( driver == NULL ) {
	  ready = false;
	  return;
	}
	  
	dataset = driver->Create(filename.c_str(), num_cols, num_rows, 
				 band_count, pixel_type,
				 options);
	// Setup georeferencing
	double geotransform[6];
	geotransform[0] = ulx;
	geotransform[3] = uly;
	geotransform[1] = 
	

	ready = true;

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
		GDALClose(dataset);
		dataset = 0;
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
	return dataset->GetRasterYSize();
}

int ProjectedRaster::getColCount()
{
	return dataset->GetRasterXSize();
}


GDALDataType ProjectedRaster::getPixelType()
{
	GDALRasterBand *band = 0;
	GDALDataType type = GDT_Unknown;
	
	band = dataset->GetRasterBand(1);
	if (band == 0) {
		return type;
	} 
	
	return band->GetRasterDataType();
}

double ProjectedRaster::getPixelSize()
{
	double geo[6];
	geo[1] = -1;
	if (isReady()) {
		dataset->GetGeoTransform(geo);
	}
		
	return geo[1];
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

	GetGDALDriverManager()->GetDriverByName("MEM");

	if (driver == 0) {
		return false;
	}

//	dataset = driver->Create(

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

