#ifndef GEOGRAPHIC_H

#define GEOGRAPHIC_H



#include "projection.h"



class Geographic : public Projection

{

public:

	Geographic();

	Geographic(double gctpParams[], ProjUnit units, ProjDatum dat);

		

private:

	//! See documentation for Projection

	void _init();



	//! See documentation for Projection

	void _forward(double lon, double lat);



	//! See documentation for Projection

	void _inverse(double x, double y);

}

