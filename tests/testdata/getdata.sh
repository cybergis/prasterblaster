#!/bin/sh

wget -c http://usgs-ybother.srv.mst.edu/static/veg_geographic_1deg.tif
wget -c http://usgs-ybother.srv.mst.edu/static/holdnorm_geographic_30min.tif
wget -c http://usgs-ybother.srv.mst.edu/static/glc_geographic_30sec.tif.bz2
bunzip2 glc_geographic_30sec.tif.bz2

wget -c http://usgs-ybother.srv.mst.edu/static/nlcd2006_landcover_4-20-11_se5.tif.bz2
bunzip2 nlcd2006_landcover_4-20-11_se5.tif.bz2


