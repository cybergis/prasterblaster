#!/bin/bash

OLDDIR=$(pwd)
DATADIR=$(dirname $0)

if command -v wget 2>/dev/null; then
    WGET="wget -c"
else
    WGET="ftp -C"
fi

cd $DATADIR

if [ -e veg_geographic_1deg.tif ]
then 
    echo "Found veg_geographic_1deg.tif"
else
    $WGET  http://usgs-ybother.srv.mst.edu/static/veg_geographic_1deg.tif
fi

if [ -e holdnorm_geographic_30min.tif ]
then 
    echo "Found holdnorm_geographic_30min.tif"
else 
    $WGET  http://usgs-ybother.srv.mst.edu/static/holdnorm_geographic_30min.tif
fi

if [ -e glc_geographic_30sec.tif ]
then
    echo "Found glc_geographic_30sec.tif"
else
    $WGET  http://usgs-ybother.srv.mst.edu/static/glc_geographic_30sec.tif.bz2
    bunzip2 glc_geographic_30sec.tif.bz2
fi

if [ -e nlcd2006_landcover_4-20-11_se5.tif ]
then
    echo "Found nlcd2006_landcover_4-20-11_se5.tif"
else 
    $WGET  http://usgs-ybother.srv.mst.edu/static/nlcd2006_landcover_4-20-11_se5.tif.bz2
    bunzip2 nlcd2006_landcover_4-20-11_se5.tif.bz2
fi

# Get control rasters for output comparison
if [ -e veg_mollweide_1deg.tif ]
then
    echo "Found veg_mollweide_1deg.tif"
else
    $WGET  http://usgs-ybother.srv.mst.edu/static/veg_mollweide_1deg.tif
fi

if [ -e holdnorm_mollweide_30min.tif ]
then
    echo "Found holdnorm_mollweide_30min.tif"
else
    $WGET  http://usgs-ybother.srv.mst.edu/static/holdnorm_mollweide_30min.tif
fi

if [ -e glc_mollweide_30sec.tif ]
then
    echo "Found glc_mollweide_30sec.tif"
else
    $WGET  http://usgs-ybother.srv.mst.edu/static/glc_mollweide_30sec.tif
fi

cd $OLDDIR
