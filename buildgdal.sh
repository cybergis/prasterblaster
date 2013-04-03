#!/bin/sh

OLDDIR=$(pwd)
PRBROOT=$(dirname $0)

cd $PRBROOT
PRBROOT=$(pwd)
cd $OLDDIR

mkdir $PRBROOT/src/gdal
cd $PRBROOT/src/gdal/

# Build proj4
wget -c http://download.osgeo.org/proj/proj-4.8.0.tar.gz
rm -rf proj-4.8.0
tar xfz proj-4.8.0.tar.gz
cd proj-4.8.0
./configure --prefix=$PRBROOT/src/gdal/
make -j2 install

# Build GDAL
cd $PRBROOT/gdal/
wget -c http://download.osgeo.org/gdal/gdal-1.9.2.tar.gz
rm -rf gdal-1.9.2
tar xfz gdal-1.9.2.tar.gz
cd gdal-1.9.2/
./configure --prefix=$PRBROOT/src/gdal/ --with-libtiff=internal --with-geotiff=internal
make -j2
make install


cd $OLDDIR
