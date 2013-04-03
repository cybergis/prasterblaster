
AM_CLFLAGS =  -W -Wall -Wcast-align -Wstrict-prototypes -Wmissing \
		-prototypes -Wpointer-arith -Wshadow -Wsign-compare \
		-Wformat=2 -Wno-format-y2k -Wimplicit -Wmissing-braces \
		-Wnested-externs -Wparentheses -Wtrigraphs -g

bin_PROGRAMS = tests systemtest
tests_CXXFLAGS  = -I../src/ -lgctp_cpp -lgtest
tests_LDFLAGS = $(all_libraries) -lgctp_cpp -lgtest -lgtest_main -L../src/ -L../src/gctp_cpp/ -L../src/gtest/ -L../src/gdal/lib/
tests_CPPFLAGS = -I ../src/gtest/include -I../src/gdal/include
tests_LDADD = -lmpi ../src/librasterblaster.la ../src/gctp_cpp/libgctp_cpp.la
tests_SOURCES = check-reprojection_tools.cc check-rastercoordtransformer.cc

systemtest_CXXFLAGS  = -I../src/ -lgctp_cpp -lgtest
systemtest_LDFLAGS = $(all_libraries) -lgctp_cpp -lgtest -L../src/ -L../src/gctp_cpp/ -L../src/gtest/ -L../src/gdal/lib/
systemtest_CPPFLAGS = -I ../src/gtest/include -I../src/gdal/include
systemtest_LDADD = -lmpi ../src/librasterblaster.la ../src/gctp_cpp/libgctp_cpp.la
systemtest_SOURCES = systemtest.cc ../src/demos/prasterblaster-pio.cc ../src/demos/sptw.cc

