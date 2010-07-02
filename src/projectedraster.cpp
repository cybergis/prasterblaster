
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

RasterChunk::~RasterChunk()
{
	delete projection;
	free(data);
	return;
	
}

ProjectedRaster::ProjectedRaster(string _filename)
{
	OGRSpatialReference sr;
	char *ref = 0;
	char **ugh = 0;
	double *params = 0;
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

		dataset = (GDALDataset*)GDALOpen( _filename.c_str(), GA_Update );

		if (dataset == 0) {
			ready = false;
			return;
		} 

		cols = dataset->GetRasterXSize();
		rows = dataset->GetRasterYSize();
		
		
		ref = strdup(dataset->GetProjectionRef());

		if (ref == 0) {
			ready = false;
			return;
		}

		double geo[6];
		geo[1] = -1;
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

		sr.exportToUSGS(&projsys, &zone, &params, &datum);
		projection = Transformer::convertProjection((ProjCode)projsys);
		if (projection == 0) {
			ready = false;
			return;
		}
		if ((ProjCode)projsys == _UTM)
			((UTM*)projection)->setZone(zone);
		projection->setParams(params);
		projection->setDatum((ProjDatum)datum);
		projection->setUnits(METER);

//		free(ref);
		ready = true;
		
	} else { // Filename specifies imagine binary raster
		bool status = loadImgRaster(filename);
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
	band_count = _band_count;
	rows = num_rows;
	cols = num_cols;
	pixel_size = _pixel_size;
	type = pixel_type;

	GDALAllRegister();

	projection = proj;

	driver = GetGDALDriverManager()->GetDriverByName(format);

	if( driver == NULL ) {
		ready = false;
		return;
	}
	  
	// Set options
	options = CSLSetNameValue( options, "INTERLEAVE", "PIXEL" );
	options = CSLSetNameValue( options, "BIGTIFF", "NO" );
	options = CSLSetNameValue( options, "TILED", "NO" );
	options = CSLSetNameValue( options, "COMPRESS", "NONE" );
	options = CSLSetNameValue( options, "PHOTOMETRIC", "MINISBLACK");

	dataset = driver->Create(filename.c_str(), num_cols, num_rows, 
				 band_count, pixel_type,
				 options);

	
	if (dataset == 0 || proj == 0) {
		ready = false;
		return;
	}

	dataset->SetGeoTransform(geotransform);
	
	band = dataset->GetRasterBand(1);


	// Setup georeferencing
	OGRSpatialReference srs;
	long zone = -1;

	/*
	if (projection->number() == _UTM) {
		zone = (int)((UTM*)projection)->zone();
	}
	srs.importFromUSGS((long)projection->number(), 
			    zone,
			    projection->params(),
			   (long)projection->datum());
	srs.Fixup();
	srs.exportToWkt(&wkt);
	printf("WKT YO: %s\n", wkt);
	dataset->SetProjection(wkt);
	CPLFree(wkt);
	*/
	
//	srs.SetProjCS(projection->name().c_str());
//	srs.SetWellKnownGeogCS( "EPSG:4052" );
	srs.importFromUSGS(proj->number(), 0, proj->params(), proj->datum());

	srs.exportToWkt(&wkt);
	dataset->SetProjection(wkt);
	CPLFree(wkt);
	ready = true;

	if (options != 0)
		CSLDestroy(options);

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
	char **options = 0; 
	GByte raster[8*8];
	char *wkt = 0;
	double geotransform[6] = { 444720, 30, 0, 3751320, 0, -30 };


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

 	FindMinBox(input, output_proj, pixel_size, ul_x, ul_y, lr_x, lr_y);
	rows = (ul_y-lr_y) / input->getPixelSize();
	cols = (lr_x-ul_x) / input->getPixelSize();

	GDALAllRegister();

	driver = GetGDALDriverManager()->GetDriverByName(format);

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

//		srs.SetProjCS(projection->name().c_str());
//		srs.SetWellKnownGeogCS( "EPSG:4052" );
		srs.importFromUSGS(output_proj->number(), 0, 
				   output_proj->params(), output_proj->datum());
		srs.SetLinearUnits("SRS_UL_METER", 1);
//		projsys = output_proj->number();
//		datum = output_proj->datum();
//		srs.importFromUSGS(projsys, zone, output_proj->params(), datum);
		srs.Fixup();
		srs.exportToWkt(&wkt);
		CPLErr err = dataset->SetProjection(wkt);
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

bool ProjectedRaster::readVector(int firstRow, int numRows, 
						   vector<unsigned char>* data)
{
	vector<unsigned char> *buf;

	buf = new vector<unsigned char>(numRows * cols * GDALGetDataTypeSize(type)/8);
	
	if(readRaster(firstRow, numRows, &(*buf)[0]) == false) {
		delete buf;
		buf = 0;
	}	
	
	return buf;

}

bool ProjectedRaster::writeVector(int firstRow, int numRows, 
						    vector<unsigned char>* data)
{
	if (writeRaster(firstRow, numRows, &(*data)[0]) == false) {
		return false;
	}


}

RasterChunk* ProjectedRaster::makeChunk(int firstRow, int numRows)
{
	RasterChunk *rc = new RasterChunk;

	rc->num_rows = numRows;
	rc->num_cols = cols;
	rc->ul_x = ul_x + pixel_size * firstRow;
	rc->ul_y = ul_y + pixel_size * firstRow;
	rc->pixel_size = pixel_size;
	rc->type = type;
	rc->projection = projection->copy();

	rc->data = (char*)malloc(band_count * rc->num_rows * rc->num_cols * (bitsPerPixel()/8));

	readRaster(firstRow, numRows, rc->data);

	if (rc->data == 0) {
		delete rc;
		delete projection;
		return 0;
	}

		
			  

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

