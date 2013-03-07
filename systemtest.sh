#!/bin/sh

./prasterblasterpio --t_srs '+proj=moll +a=6370997 +b=6370997' -r mean tests/testdata/glc_geographic_30sec.tif tests/testdata/glc_mollweide_30sec.tif

./prasterblasterpio --t_srs '+lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m' -r mean \
                  'tests/testdata/nlcd2006_landcover_4-20-11_se5.tif' 'tests/testoutput/nlcd2006_landcover_4-20-11_se5_mollweide.tif'


