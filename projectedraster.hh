

#ifndef RASTER_HH
#define RASTER_HH

#include <string>

#include "gctp_cpp/constants.h"
#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"

#include "rasterinfo.h"

using namespace std; // Don't do this :(

class ProjectedRaster
{
public:
	ProjectedRaster();
	ProjectedRaster(string filename);
	ProjectedRaster(long num_rows, long num_cols, long pixel_bits);
	~ProjectedRaster();
	void* getData();
	Projection* getProjection();
        bool isReady();
	ProjectedRaster* getSubraster(int firstrow, int lastrow);
	
	// Area
	void setUL(double ul_x, double ul_y);
	void setRowCount(int rows);
	void setColCount(int cols);
	int getRowCount();
	int getColCount();
	
	// Pixel description
	void setDataType(const std::string &datatype);
	void setPixelSize(double pixsize);
	void setBitCount(int bits);
	void setSigned();
	void setUnsigned();
	void setInteger();
	void setFloat();
	string getDataType();
	double getPixelSize();
	int getBitCount();
	bool isSigned();
	bool isInteger();
	bool isFloat();
	
	
	// Projection
	void setProjection(ProjCode p);
	void setZoneNumber(int zone);
	void setUnit(ProjUnit unit);
	void setDatum(ProjDatum pd);
	void setGctpParams(double params[]);
	int getZoneNumber();
	ProjDatum getDatum();
	double* getGctpParams();

	// Reprojection
	void setReprojection();

	// Pixel Values
	void getValueGeo(double lat, double lon, void* val);
	

	// Members
	double ul_x, ul_y;

	// Pixel Description
	int rows, cols, bitCount;
	bool issigned;
	bool integer;
	double pixsize; // In degrees
	std::string type;

	// Projection
	int zoneNumber;
	ProjCode projectionCode;
	ProjDatum datum;
	ProjUnit unit;
	double gctpParams[16];
	
	void* data;
private:
	bool readImgRaster(string filename);
	bool readRaster(string filename);
	Projection *projection;
	Transformer t;
        ProjectedRaster& operator=(ProjectedRaster& a);
	bool ready;
};

#endif //RASTER_HH
