


noinst_LTLIBRARIES = librasterblaster.la
librasterblaster_la_SOURCES = projectedraster.cc \
                             configuration.cc \
                             quadtree.cc \
                             rastercoordtransformer.cc \
                             reprojection_tools.cc
#librasterblaster_a_LDFLAGS = $(all_libraries)  -Lgctp_cpp
librasterblaster_la_LIBADD =  gctp_cpp/libgctp_cpp.la
librasterblaster_la_LDFLAGS = -lgdal

AM_CPPFLAGS = $(all_includes) -IpRPL -Igctp_cpp

AM_CXXFLAGS = -W -Wall -Wextra -Wcast-align \
		-Wpointer-arith -Wsign-compare -g\
		-Wformat=2 -Wno-format-y2k  -Wmissing-braces -Wparentheses \
		-Wtrigraphs -Wstrict-aliasing=2

SUBDIRS =  gctp_cpp . gtest demos
