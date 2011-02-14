
AM_CLFLAGS =  -W -Wall -Wcast-align -Wstrict-prototypes -Wmissing \
		-prototypes -Wpointer-arith -Wshadow -Wsign-compare \
		-Wformat=2 -Wno-format-y2k -Wimplicit -Wmissing-braces \
		-Wnested-externs -Wparentheses -Wtrigraphs -g

bin_PROGRAMS = tests
tests_CXXFLAGS  = -I../src/
tests_LDFLAGS = $(all_libraries) -L../src/
tests_LDADD = -lmpi -lrasterblaster -lgtest
tests_SOURCES = check_reprojector.cpp

