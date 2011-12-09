SUBDIRS = gctp_cpp gtest


noinst_LIBRARIES = librasterblaster.a
librasterblaster_a_SOURCES = projectedraster.cpp \
                 rasterinfo.cpp \
	         reprojector.cpp \
                 tinystr.cpp \
                 tinyxml.cpp \
                 tinyxmlerror.cpp \
                 tinyxmlparser.cpp \
                 rasterxml.cpp \
		 resampler.cpp \
		 driver.cpp
#librasterblaster_a_LDFLAGS = $(all_libraries)  -Lgctp_cpp
librasterblaster_a_LIBADD =  gctp_cpp/libgctp_cpp.a

bin_PROGRAMS = prasterblaster
prasterblaster_CXXFLAGS = -lgdal -g
prasterblaster_LDFLAGS = $(all_libraries)  -Lgctp_cpp
prasterblaster_LDADD = librasterblaster.a gctp_cpp/libgctp_cpp.a 
prasterblaster_SOURCES =  main.cpp 

AM_CPPFLAGS = $(all_includes) -IpRPL -Igctp_cpp

AM_CXXFLAGS = -W -Wall -Wextra -Wcast-align \
		-Wpointer-arith -Wsign-compare -g\
		-Wformat=2 -Wno-format-y2k  -Wmissing-braces -Wparentheses \
		-Wtrigraphs -std=c++0x # -Wshadow -Wstrict-prototypes -Wnested-externs -Wimplicit
