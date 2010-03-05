SUBDIRS = pRPL gctp_cpp

bin_PROGRAMS = prasterblaster
prasterblaster_LDFLAGS = $(all_libraries) -LpRPL -Lgctp_cpp
prasterblaster_LDADD = -lprp -lgctp_cpp
prasterblaster_SOURCES =  projectedraster.cpp \
                 rasterinfo.cpp \
	         reprojector.cpp \
                 tinystr.cpp \
                 tinyxml.cpp \
                 tinyxmlerror.cpp \
                 tinyxmlparser.cpp \
                 rasterxml.cpp \
                 driver.cpp \
                 resampler.cpp

AM_CPPFLAGS = $(all_includes) -IpRPL -Igctp_cpp
