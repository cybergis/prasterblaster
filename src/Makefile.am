SUBDIRS = gctp_cpp gtest


noinst_LIBRARIES = librasterblaster.a
librasterblaster_a_SOURCES = projectedraster.cc \
                             configuration.cc \
                             projectedraster.cc \
                             quadtree.cc \
                             rastercoordtransformer.cc \
                             reprojection_tools.cc
#librasterblaster_a_LDFLAGS = $(all_libraries)  -Lgctp_cpp
librasterblaster_a_LIBADD =  gctp_cpp/libgctp_cpp.a

AM_CPPFLAGS = $(all_includes) -IpRPL -Igctp_cpp

AM_CXXFLAGS = -W -Wall -Wextra -Wcast-align \
		-Wpointer-arith -Wsign-compare -g\
		-Wformat=2 -Wno-format-y2k  -Wmissing-braces -Wparentheses \
		-Wtrigraphs -Wstrict-aliasing=2
