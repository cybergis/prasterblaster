# Filename: Makefile.am
# Author: David Mattli <dmattli@usgs.gov>
# License: PUBLIC DOMAIN

noinst_LTLIBRARIES = libgctp_cpp.la

AM_CXXFLAGS = -W -Wall -Wextra -Wcast-align \
                -Wpointer-arith -Wsign-compare -g\
                -Wformat=2 -Wno-format-y2k  -Wmissing-braces -Wparentheses \
                -Wtrigraphs -Wstrict-aliasing=2

libgctp_cpp_la_CXXFLAGS = -fPIC -L../gdal/lib
libgctp_cpp_la_CPPFLAGS = -I../gdal/include/
libgctp_cpp_la_SOURCES = alaskaconformal.cpp \
           albersConEqArea.cpp \
           azequidistant.cpp \
           constants.cpp \
           coordinate.cpp \
           cproj.cpp \
           equidistantc.cpp \
           equirectangular.cpp \
           genvertnsp.cpp \
           gnomonic.cpp \
           goodeh.cpp \
           hammer.cpp \
           hotinobmerc.cpp \
           intmollweide.cpp \
           lambertazimuth.cpp \
           lambertcc.cpp \
           mercator.cpp \
           miller.cpp \
           Mollweide.cpp \
           oblatedeqarea.cpp \
           orthographic.cpp \
           polarstereo.cpp \
           polyconic.cpp \
           projection.cpp \
           robinson.cpp \
           sinusoidal.cpp \
           spaceobmerc.cpp \
           stereo.cpp \
           transformer.cpp \
           transversemercator.cpp \
           util.cpp \
           utm.cpp \
           vandergrinten.cpp \
           wagneriv.cpp \
           wagnervii.cpp 
