
#include <QString>

#include <string>
#include <cstdio>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "rasterreader.hh"
#include "rasterinfo.h"
#include "projectedraster.hh"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/projection.h"
#include "gctp_cpp/utm.h"

using namespace std;

ProjectedRaster* RasterReader::readRaster(std::string filename)
{
  GDALDataset  *dataset;
  GDALRasterBand *rasterband;
  int pixelsize, xsize, ysize, bandcount;
  long projcode, projzone, projdatum;
  double *params = (double*)calloc(16, sizeof(double));
  ProjectedRaster* raster = new ProjectedRaster;
  Projection *proj = 0;
  Transformer t;
  char* temp;
  void* data;

  GDALAllRegister();
  
  dataset = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly );
  if( dataset == NULL ){
    return 0;
  }
  
  rasterband = dataset->GetRasterBand(1);
  pixelsize =  GDALGetDataTypeSize(rasterband->GetRasterDataType())/8;
  xsize = dataset->GetRasterXSize();
  ysize = dataset->GetRasterYSize();
  bandcount = dataset->GetRasterCount();
  OGRSpatialReference spatialref(dataset->GetProjectionRef());
  spatialref.Fixup();
  spatialref.exportToPrettyWkt(&temp);
  printf("WKT!: %s\n", temp);
  spatialref.exportToUSGS(&projcode, &projzone, &params, &projdatum);
  printf("ProjCode: %ld, ProjZone: %ld, ProjDatum: %ld\n", projcode,
	 projzone, projdatum);
  // Setup projection info
  proj = t.convertProjection((ProjCode)projcode);
  if (proj != 0) {
    proj->setDatum((ProjDatum)projdatum);
    proj->setParams(params);
  } else {
    printf("invalid projection!\n");
  }

  // Read Raster Data
  data = malloc(pixelsize*xsize*ysize*bandcount);
  if (data != 0) {
    // TODO: Instead of NULL, use a map to ensure ARGB order
    CPLErr err = dataset->RasterIO(GF_Read, 0, 0, xsize, ysize, data,
		      xsize, ysize, rasterband->GetRasterDataType(),
		      bandcount, NULL, bandcount, xsize*bandcount, 1);
    if (err != CE_None) {
      printf("Error reading image!\n");
    } else {
      printf("Success reading image!\n");
    }
  } else {
    printf("Error allocating memory!\n");
  }
  
  raster->data = data;

  GDALClose(dataset);

  return raster;
}

ProjectedRaster* RasterReader::readImgRaster(std::string filename, int rank, int num_procs)
{
  string imgname, xmlname;
  imgname = xmlname = filename;
  imgname.append(".img");
  xmlname.append(".xml");
  RasterInfo in_info(xmlname.c_str());
  int fd = -1;
  int readsize = 0;
  int rastersize = in_info.rows()/num_procs * in_info.cols();
  int overflow = (in_info.rows() % num_procs) * in_info.cols();
  ProjectedRaster *pr = new ProjectedRaster;

  if (rank != num_procs-1) {
    overflow = 0;
  }

  pr->data = calloc(rastersize + overflow, in_info.bitCount()/8);
  if (pr->data == 0) {
    printf("Data allocation failed\n");
    return 0;
  }
  // Read in image data
  errno = 0;
  fd = open(imgname.c_str(), O_RDONLY);
  lseek(fd, rank * rastersize, SEEK_SET);
  readsize = read(fd, pr->data, 
		  (rastersize+overflow)*(in_info.bitCount()/8));
  if (readsize != (rastersize * (in_info.bitCount()/8))) {
    printf("Read error: %d!\n", readsize);
    printf("Supposed to be: %d\n", rastersize*(in_info.bitCount()/8));
    printf("ERROR: %s\n", strerror(errno));
    return 0;
  }
  close(fd);

  // Setup projection
  int num_rows = in_info.rows()/num_procs;
  if (rank == num_procs-1)
    num_rows += in_info.rows() % num_procs;
  pr->setProjection((ProjCode)in_info.projectionNumber());
  pr->setUL(in_info.ul_X(), in_info.ul_Y() 
	    - (rank * (in_info.rows()/num_procs) * in_info.pixelSize()));
  pr->setDatum((ProjDatum)in_info.datumNumber());
  pr->setRowCount(num_rows);
  pr->setColCount(in_info.cols());
  pr->setPixelSize(in_info.pixelSize());
  if (in_info.isSigned())
    pr->setSigned();
  else
    pr->setUnsigned();
  pr->setUnit((ProjUnit)in_info.unitNumber());
  pr->setGctpParams(in_info.allGctpParams());
  
  return pr;
}

void RasterReader::writeRaster(std::string filename,
		 ProjectedRaster *raster)
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
  proj = raster->getProjection();
  if (raster->getProjection()->number() == _UTM) {
    sr.importFromUSGS((long)proj->number(), ((UTM*)proj)->zone(),
		      proj->params(), (long)proj->datum());
  } else {
    sr.importFromUSGS((long)proj->number(),  0,
		      proj->params(), (long)proj->datum());
  }
  sr.SetLinearUnits(SRS_UL_METER, 1);
  sr.exportToWkt( &wkt );

  // Create Dataset
  //  options = CSLSetNameValue(options, "INTERLEAVE", "PIXEL");
  //  options = CSLSetNameValue(options, "COMPRESS", "PACKBITS");
  ds = driver->Create(filename.c_str(), raster->getColCount(),
    		      raster->getRowCount(), raster->getBitCount()/8, 
		      GDT_Byte, options);
  


  if (ds != 0) {
    // TODO: Support more than one band
    printf("Writing raster: %d cols %d rows\n", raster->getColCount(), raster->getRowCount());
    if (ds->RasterIO(GF_Write, 0, 0, raster->getColCount(), raster->getRowCount(),
		     raster->getData(), raster->getColCount(), raster->getRowCount(),
		     GDT_Byte, 1, 0, 0, 0, 0) == CE_Failure) {
      printf("RasterIO failed...\n");
    }

    ds->SetProjection(wkt);
    GDALClose(ds);
  }
  return;
}

