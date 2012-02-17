
AM_CLFLAGS =  -W -Wall -Wcast-align -Wstrict-prototypes -Wmissing \
		-prototypes -Wpointer-arith -Wshadow -Wsign-compare \
		-Wformat=2 -Wno-format-y2k -Wimplicit -Wmissing-braces \
		-Wnested-externs -Wparentheses -Wtrigraphs -g

bin_PROGRAMS = tests
tests_CXXFLAGS  = -I../src/ -lgctp_cpp -lgtest -std=c++0x
tests_LDFLAGS = $(all_libraries) -lgctp_cpp -lgtest -L../src/ -L../src/gctp_cpp/ -L../src/gtest/
tests_CPPFLAGS = -I ../src/gtest/include
tests_LDADD = -lmpi ../src/librasterblaster.a ../src/gctp_cpp/libgctp_cpp.a -lgdal
tests_SOURCES = check_driver.cpp check_reprojector.cpp check_projectedraster.cpp #check_chunker.cpp \
		check_rastercache.cpp raster_tools.cpp 

