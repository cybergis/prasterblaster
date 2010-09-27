SUBDIRS = gctp_cpp 

AM_CFLAGS = -W -Wall -Wcast-align -Wstrict-prototypes -Wmissing \
		-prototypes -Wpointer-arith -Wshadow -Wsign-compare \
		-Wformat=2 -Wno-format-y2k -Wimplicit -Wmissing-braces \
		-Wnested-externs -Wparentheses -Wtrigraphs

lib_LTLIBRARIES = librasterblaster.la
librasterblaster_la_SOURCES = projectedraster.cpp \
                 rasterinfo.cpp \
	         reprojector.cpp \
                 tinystr.cpp \
                 tinyxml.cpp \
                 tinyxmlerror.cpp \
                 tinyxmlparser.cpp \
                 rasterxml.cpp \
		 resampler.cpp
librasterblaster_la_LDFLAGS = $(all_libraries)  -Lgctp_cpp
librasterblaster_la_LIBADD =  -lgctp_cpp -lmpi

bin_PROGRAMS = prasterblaster
prasterblaster_LDFLAGS = $(all_libraries)  -Lgctp_cpp
prasterblaster_LDADD =  -lgctp_cpp -lmpi -lrasterblaster
prasterblaster_SOURCES =  driver.cpp

AM_CPPFLAGS = $(all_includes) -IpRPL -Igctp_cpp
