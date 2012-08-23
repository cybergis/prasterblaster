/*!
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE 
 *
 * This software is in the public domain, furnished "as is", without
 * technical support, and with no warranty, express or implied, as to
 * its usefulness for any purpose.
 *
 * @section DESCRIPTION
 *
 * The ProjectedRaster class represents a raster with a location and a projection.
 *
 */

#include <string>
#include <sstream>
#include <cmath>
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

#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mercator.h"
#include "gctp_cpp/constants.h"
#include "gctp_cpp/utm.h"
#include "reprojector.hh"
#include "rasterchunk.hh"
#include "sharedptr.hh"

#include "projectedraster.hh"

using std::string;

ProjectedRaster::ProjectedRaster(string _filename)
{
	OGRSpatialReference sr;

	filename = _filename;
	dataset = 0;
	data = 0;
	rows = cols = -1;

	GDALAllRegister();

	ready = loadRaster(filename);

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



bool ProjectedRaster::CreateRaster(string _filename,
				   shared_ptr<ProjectedRaster> input,
				   shared_ptr<Projection> output_proj,
				   GDALDataType pixel_type,
				   double _pixel_size)
{
	bool status = false;
	double ulx, uly;
	int num_cols, num_rows;

	Area out_area =  ProjectedMinbox(input,
					 output_proj);
				
	ulx = out_area.ul.x;
	uly = out_area.ul.y;

	num_rows = (int)(ceil(out_area.ul.y - out_area.lr.y) / _pixel_size);
	num_cols = (int)(ceil(out_area.lr.x - out_area.ul.x) / _pixel_size);

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

	if (isReady() && dataset != 0) { // GTiff
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
	
	
	// TODO: Verify parameters!
	if (isReady() && dataset != 0) { // GTiff
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
}

RasterChunk::RasterChunk* ProjectedRaster::createRasterChunk(Area area)
{
	RasterChunk::RasterChunk *temp = createAllocatedRasterChunk(area);


	if (temp == NULL) {
		return NULL;
	}

        // Read area of raster
        if (isReady() && dataset != 0) {
                if (dataset->RasterIO(GF_Read,
                                      area.ul.x,
                                      area.ul.y,
				      temp->column_count_,
				      temp->row_count_,
				      temp->pixels_,
				      temp->column_count_,
				      temp->row_count_,
                                      temp->pixel_type_,
                                      bandCount(),
                                      NULL,
                                      0,0,0) != CE_None)
                {
                        // Cleanup and return error
                        delete temp;
                        return NULL;
                }
                                      

        } else{ 
		delete temp;
		return NULL;
	}

	return temp;
}

RasterChunk::RasterChunk* ProjectedRaster::createAllocatedRasterChunk(Area area)
{
	RasterChunk::RasterChunk *temp = createEmptyRasterChunk(area);	
	size_t buffer_size = (temp->row_count_ * temp->column_count_);
	unsigned char *pixels = (unsigned char*)calloc(buffer_size, GDALGetDataTypeSize(getPixelType())/8);

	if (pixels == NULL) {
		fprintf(stderr, "Allocation error!\n");
		delete temp;
		return NULL;
	}

	temp->pixels_ = pixels;

	return temp;
	

}

RasterChunk::RasterChunk* ProjectedRaster::createEmptyRasterChunk(Area area)
{
	RasterChunk::RasterChunk *temp = new RasterChunk::RasterChunk;

	temp->projection_ = shared_ptr<Projection>(getProjection());
	temp->raster_location_ = area.ul;
	temp->ul_projected_corner_ = Coordinate(ul_x+(area.ul.x*getPixelSize()), ul_y-(area.ul.y*getPixelSize()), UNDEF);
	temp->pixel_size_ = getPixelSize();
	temp->row_count_ = area.lr.y - area.ul.y + 1;
	temp->column_count_ = area.lr.x - area.ul.x + 1;
	temp->pixel_type_ = getPixelType();
	temp->band_count_ = band_count;
	temp->pixels_ = NULL;

	return temp;
}

bool ProjectedRaster::writeRasterChunk(RasterChunk::RasterChunk *chunk)
{
	// TODO: Add some more checks

	if (dataset->RasterIO(GF_Write,
			      chunk->raster_location_.x,
			      chunk->raster_location_.y,
			      chunk->column_count_,
			      chunk->row_count_,
			      chunk->pixels_,
			      chunk->column_count_,
			      chunk->row_count_,
			      chunk->pixel_type_,
			      chunk->band_count_,
			      NULL, 0, 0, 0) != CE_None) {
		// Error!
		fprintf(stderr, "Error writing RasterChunk %p\n", chunk->pixels_);
		return false;

	}
	dataset->FlushCache();
			      

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
		fprintf(stderr, "Error opening GDAL dataset\n");
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
	projection->setParams(p);
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
	options = CSLSetNameValue( options, "BIGTIFF", "IF_SAFER" );
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
	geotransform[1] = _pixel_size;
	geotransform[5] = -_pixel_size;
	_dataset->SetGeoTransform(geotransform);

	srs.importFromUSGS(_projection->number(), 0, _projection->params(), _projection->datum());

	char *wkt = NULL;
	srs.exportToWkt(&wkt);
	if (wkt != NULL) {
		_dataset->SetProjection(wkt);
		OGRFree(wkt);
	} else {
		fprintf(stderr, "\nERROR: Unsupported projection: %s\n", _projection->name().c_str());
	}

	GDALClose(_dataset);

	if (options != 0)
		CSLDestroy(options);

	return true;

}
