
/* 
 * Filename: check_reprojector.cpp
 * Author: David Mattli <dmattli@usgs.gov>
 * License: PUBLIC DOMAIN
 */


#include <algorithm>
#include <cfloat>
#include <string>
#include <vector>


#include <boost/shared_ptr.hpp>

#include <gtest/gtest.h>
#include <mpi.h>

#include "gctp_cpp/projection.h"
#include "gctp_cpp/transformer.h"
#include "gctp_cpp/mollweide.h"

#include "driver.hh"
#include "reprojector.hh"
#include "prog_state.hh"


using std::string;
using boost::shared_ptr;


class MinboxTest : public ::testing::Test {
protected:
	virtual void SetUp() {

	}

	

};

void updateMinbox(double x, double y, Area *minbox) 
{
	if (x > -DBL_MAX && x < minbox->ul.x) 
		minbox->ul.x = x;
	if (x < DBL_MAX && x > minbox->lr.x) 
		minbox->lr.x = x;
	if (y > -DBL_MAX && y < minbox->lr.y) 
		minbox->lr.y = y;
	if (y < DBL_MAX && y > minbox->ul.y) 
		minbox->ul.y = y;
	
	return;
}


TEST(MinboxTest, geoextent) {
	Projection *p = 0;
	Area minbox;
	Area geominbox;
	double temp1, temp2; 	
	geominbox.ul.x = minbox.ul.x = DBL_MAX;
	geominbox.lr.x = minbox.lr.x = -DBL_MAX;
	geominbox.ul.y = minbox.ul.y = -DBL_MAX;
	geominbox.lr.y = minbox.lr.y = DBL_MAX;

	p = Transformer::convertProjection(POLYC);
	
	ASSERT_NE(p, (Projection*)0);
	
	p->setUnits(METER);
	p->setDatum(WGS_84);
	
	for (double lon = -180.0; lon <= 180.0; lon += 0.5) {
		for (double lat = -90; lat <= 90; lat += 0.5) {
			p->forward(lon, lat);
			ASSERT_FALSE(p->errorOccured());
			updateMinbox(p->x(), p->y(), &minbox);
		}
	}

	for (double x = minbox.ul.x; x <= minbox.lr.x; x += 10000) {
		for (double y = minbox.lr.y; y <= minbox.ul.y; y += 10000) {
			p->inverse(x, y);
			ASSERT_FALSE(p->errorOccured());
			updateMinbox(p->lon(), p->lat(), &geominbox);
		}
		
	}
	

	p->inverse(minbox.ul.x, minbox.ul.y, &geominbox.ul.x, &geominbox.ul.y);
	p->inverse(minbox.lr.x, minbox.lr.y, &geominbox.lr.x, &geominbox.lr.y);

	printf("Calculated Minbox: Proj:%f %f %f %f\n"
	       "Geo: %f %f %f %f\n",
	       minbox.ul.x, minbox.ul.y, minbox.lr.x, 
	       minbox.lr.y,
	       geominbox.ul.x, geominbox.ul.y, geominbox.lr.x, geominbox.lr.y);
}
