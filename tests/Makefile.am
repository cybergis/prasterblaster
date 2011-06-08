
AM_CLFLAGS =  -W -Wall -Wcast-align -Wstrict-prototypes -Wmissing \
		-prototypes -Wpointer-arith -Wshadow -Wsign-compare \
		-Wformat=2 -Wno-format-y2k -Wimplicit -Wmissing-braces \
		-Wnested-externs -Wparentheses -Wtrigraphs -g

bin_PROGRAMS = tests
tests_CXXFLAGS  = -I../src/
tests_LDFLAGS = $(all_libraries) -L../src/ -L../src/gtest/
tests_CPPFLAGS = -I ../src/gtest/include
tests_LDADD = -lmpi -lrasterblaster -lgtest
tests_SOURCES = check_reprojector.cpp check_projectedraster.cpp check_chunker.cpp \
		check_rastercache.cpp

