
#include <QString>

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "rasterinfo.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mercator.h"
#include "reprojector.hh"

#include "projectedraster.hh"

using namespace std;

ProjectedRaster::ProjectedRaster()
{
  data = 0;
  ready = false;
  ul_x = ul_y = 0;
  rows = cols = 0;
  bitCount = 8;
  issigned = false;
  integer = true;
  projection = 0;

  for(int i = 0; i < 16; ++i) {
    gctpParams[i] = 0;
  }

  return;
}

ProjectedRaster::ProjectedRaster(long num_rows, long num_cols, 
				 long pixel_bits)
{
  /*
  data = 0;
  setProjection(proj);
  if (input.isSigned())
    setSigned();
  else
    setUnsigned();
  setBitCount(input.getBitCount());
  if(input.isInteger())
    setInteger();
  else 
    setFloat();
  setDatum(input.getDatum());
  setPixelSize(input.getPixelSize());
  
  
  // Calculate minbox and allocate memory
  FindMinBox(input, getProjection(), getPixelSize(),
	     ul_x, ul_y, lr_x, lr_y);

  cols = (lr_x - ul_x) / getPixelSize();
  rows = (ul_y - lr_y) / getPixelSize();
  */
  data = calloc(num_cols*num_rows, (pixel_bits/8));
  ready = false;
  ul_x = ul_y = 0;
  rows = num_rows;
  cols = num_cols;
  bitCount = pixel_bits;
  issigned = false;
  integer = true;
  projection = 0;

  for(int i = 0; i < 16; ++i) {
    gctpParams[i] = 0;
  }
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

void getValueGeo(double lat, double lon, void* val)
{
  
  return;
}
