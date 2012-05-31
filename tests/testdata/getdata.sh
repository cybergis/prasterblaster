#!/bin/sh

wget -c http://mattli.us/veg_mollweide_1deg.tif
wget -c http://mattli.us/holdnorm_geographic_30min.tif
wget -c http://mattli.us/glc_geographic_30sec.tif.bz2
bunzip2 glc_geographic_30sec.tif.bz2

wget -c http://mattli.us/nlcd2006_landcover_4-20-11_se5.tif.bz2
bunzip2 nlcd2006_landcover_4-20-11_se5.tif.bz2


