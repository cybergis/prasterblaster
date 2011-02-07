SUBDIRS = gctp_cpp 


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
librasterblaster_la_LIBADD =  gctp_cpp/libgctp_cpp.a

bin_PROGRAMS = prasterblaster
prasterblaster_CXXFLAGS = -lgdal
prasterblaster_LDFLAGS = $(all_libraries)  -Lgctp_cpp
prasterblaster_LDADD = librasterblaster.la 
prasterblaster_SOURCES =  driver.cpp 

AM_CPPFLAGS = $(all_includes) -IpRPL -Igctp_cpp

AM_CXXFLAGS = -W -Wall -Wcast-align -Wstrict-prototypes  \
		-Wpointer-arith -Wsign-compare \
		-Wformat=2 -Wno-format-y2k -Wimplicit -Wmissing-braces \
		-Wnested-externs -Wparentheses -Wtrigraphs # -Wshadow
