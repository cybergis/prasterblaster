


libRasterBlaster Overview
===============

Most of the functionality in librasterblaster is concerned with the
manipulation of three distinct coordinate spaces that I call: raster,
projected, and geographic.

The raster space is used to address the individual pixel values in the
raster. It is a 2D cartesian coordinate plane with (0,0) in the
top-left corner. The y-axis grows down and the x-axis grows to the
right. 

The projected coordinate space is a cartesian plane determined by the
map projection used. In projected coordinates (0,0) is in the center
of the space and the y-axis grows up.

Converting from projected coordinates to raster coordinates and back
is a simple transformation.

Geographic coordinates are the spheroidal longitude and latitude used
to describe locations on the earth. Converting projected coordinates
to geographic requires the use of map projection equations.

Coordinate spaces
==============

The projected coordinate system a raster uses is specified by its
projection string and its location in the projected space is specified
by its upper-left corner.

For example:

$ gdalinfo tests/testdata/nlcd2006_landcover_4-20-11_se5.tif |head -30

produces the output:

=================================
Driver: GTiff/GeoTIFF
Files: tests/testdata/nlcd2006_landcover_4-20-11_se5.tif
Size is 161190, 104424
Coordinate System is:
PROJCS["Albers Conical Equal Area",
    GEOGCS["NAD83",
        DATUM["North_American_Datum_1983",
            SPHEROID["GRS 1980",6378137,298.2572221010002,
                AUTHORITY["EPSG","7019"]],
            AUTHORITY["EPSG","6269"]],
        PRIMEM["Greenwich",0],
        UNIT["degree",0.0174532925199433],
        AUTHORITY["EPSG","4269"]],
    PROJECTION["Albers_Conic_Equal_Area"],
    PARAMETER["standard_parallel_1",29.5],
    PARAMETER["standard_parallel_2",45.5],
    PARAMETER["latitude_of_center",23],
    PARAMETER["longitude_of_center",-96],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]]]
Origin = (-2493045.000000000000000,3310005.000000000000000)
Pixel Size = (30.000000000000000,-30.000000000000000)
Metadata:
  AREA_OR_POINT=Area
Image Structure Metadata:
  INTERLEAVE=BAND
Corner Coordinates:
Upper Left  (-2493045.000, 3310005.000) (130d13'58.18"W, 48d42'26.63"N)
===================================

The coordinate system is specified in the WKT format and the
upper-left corner of the raster is listed on the last line.

Using this information we could convert any raster coordinate to its
equivalent projection coordinate. For example:


Where is the raster coordinate (33, 150) in projected coordinates?

X_proj = (X * pixel_size) + upper-left-x
       = (33 * 30) + -2493045 
       = -2492055

Y_proj = upper-left-y - (Y * pixel_size)
       = 3310005 - (150 * 30)
       = 326505

So the equivalent projected coordinate is (-2492055, 326505).

Finding the geographic coordinate requires a more complicated transformation, so we'll use the invproj program:

$ invproj +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs

-249205  326505
98d27'56.789"W  25d57'33.691"N <-- This is the geographic coordinate


Reprojection Process
====================

The first task is to determine the extent of the projected output
space, this combined with the size of the pixels tells us the size of
the output raster in raster coordinates.

To find the output projected space we perform a minbox operation. 

The output raster coordinate space is then partitioned. This is just a
logical partitioning, no pixel data has been read yet.

For each output raster partition we now perform another minbox
operation, to find the corresponding partition in the input raster
space.

For each partition pair we read the input, resample to the output
partition, and then write the output partition. This last operation is
performed in parallel on multiple processors. 

