


lib_LTLIBRARIES = librasterblaster.la
librasterblaster_la_SOURCES = projectedraster.cc \
                             configuration.cc \
                             quadtree.cc \
                             rastercoordtransformer.cc \
                             reprojection_tools.cc \
                             rasterchunk.cc
librasterblaster_la_LDFLAGS = $(all_libraries) -lgdal -Lgctp_cpp -Lgdal/lib/
librasterblaster_la_LIBADD =  gctp_cpp/libgctp_cpp.la

AM_CPPFLAGS = $(all_includes) -IpRPL -Igctp_cpp -Igdal/include

AM_CXXFLAGS = -W -Wall -Wextra -Wcast-align \
		-Wpointer-arith -Wsign-compare -g\
		-Wformat=2 -Wno-format-y2k  -Wmissing-braces -Wparentheses \
		-Wtrigraphs -Wstrict-aliasing=2

SUBDIRS =  gctp_cpp . gtest demos
