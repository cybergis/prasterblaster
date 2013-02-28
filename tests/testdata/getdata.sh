#!/bin/bash

OLDDIR= $(pwd)
DATADIR= $(dirname $0)

cd $DATADIR

if [ -e veg_geographic_1deg.tif ]
then 
    echo "Found veg_geographic_1deg.tif"
else
    wget -c http://usgs-ybother.srv.mst.edu/static/veg_geographic_1deg.tif
fi

if [ -e holdnorm_geographic_30min.tif ]
then 
    echo "Found holdnorm_geographic_30min.tif"
else 
    wget -c http://usgs-ybother.srv.mst.edu/static/holdnorm_geographic_30min.tif
fi

if [ -e glc_geographic_30sec.tif ]
then
    echo "Found glc_geographic_30sec.tif"
else
    wget -c http://usgs-ybother.srv.mst.edu/static/glc_geographic_30sec.tif.bz2
    bunzip2 glc_geographic_30sec.tif.bz2
fi

if [ -e nlcd2006_landcover_4-20-11_se5.tif ]
then
    echo "Found nlcd2006_landcover_4-20-11_se5.tif"
else 
    wget -c http://usgs-ybother.srv.mst.edu/static/nlcd2006_landcover_4-20-11_se5.tif.bz2
    bunzip2 nlcd2006_landcover_4-20-11_se5.tif.bz2
fi

cd $OLDDIR