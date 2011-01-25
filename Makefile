SRCDIR = src/
PRBDIR = src/
GCTPDIR = src/gctp_cpp/
OBJDIR = src/objs
TESTDIR = tests/

LCFLAGS = -fPIC $(CFLAGS) -g -Wall -W -O2 -L$(GCTPDIR) -L.
LINCLUDES = $(INCLUDES) -Isrc -I$(GCTPDIR)
LLIBS = $(LIBS) -lmpi -lgdal1.6.0 -lgctp_cpp

CXX = mpic++

_GCTPSOURCES = alaskaconformal.cpp \
           albersConEqArea.cpp \
           azequidistant.cpp \
           constants.cpp \
           coordinate.cpp \
           cproj.cpp \
           equidistantc.cpp \
           equirectangular.cpp \
           genvertnsp.cpp \
           gnomonic.cpp \
           goodeh.cpp \
           hammer.cpp \
           hotinobmerc.cpp \
           intmollweide.cpp \
           lambertazimuth.cpp \
           lambertcc.cpp \
           mercator.cpp \
           miller.cpp \
           Mollweide.cpp \
           oblatedeqarea.cpp \
           orthographic.cpp \
           polarstereo.cpp \
           polyconic.cpp \
           projection.cpp \
           robinson.cpp \
           sinusoidal.cpp \
           spaceobmerc.cpp \
           stereo.cpp \
           transformer.cpp \
           transversemercator.cpp \
           util.cpp \
           utm.cpp \
           vandergrinten.cpp \
           wagneriv.cpp \
           wagnervii.cpp

_PRBSOURCES =  projectedraster.cpp \
               rasterinfo.cpp \
               reprojector.cpp \
               tinystr.cpp \
               tinyxml.cpp \
               tinyxmlerror.cpp \
               tinyxmlparser.cpp \
               rasterxml.cpp \
               resampler.cpp

_TESTSOURCES = check_projectedraster.cpp

GCTPSOURCES = $(patsubst %,$(GCTPDIR)/%,$(_GCTPSOURCES))
PRBSOURCES = $(patsubst %,$(PRBDIR)/%,$(_PRBSOURCES))
TESTSOURCES = $(patsubst %,$(TESTDIR)/%,$(_TESTSOURCES))
GCTPOBJS = $(patsubst %.cpp,%.o,$(GCTPSOURCES))
PRBOBJS = $(patsubst %.cpp,%.o,$(PRBSOURCES))
TESTOBJS = $(patsubst %.cpp,%,$(TESTSOURCES))

all: prasterblaster libprasterblaster.so

.cpp.o:
	$(CXX) $(LCFLAGS) $(LINCLUDES) -c $< -o $@

$(GCTPDIR)/libgctp_cpp.a: $(GCTPOBJS)
	ar rcs $(GCTPDIR)/libgctp_cpp.a $(GCTPOBJS)

libprasterblaster.a: $(PRBOBJS) $(GCTPDIR)/libgctp_cpp.a
	ar rcs libprasterblaster.a $(PRBOBJS) $(GCTPOBJS)

libprasterblaster.so: $(PRBOBJS) $(GCTPDIR)/libgctp_cpp.a
	$(CXX) -shared $(LCFLAGS) $(LINCLUDES) $(PRBOBJS) $(GCTPOBJS) -o libprasterblaster.so $(LLIBS)

prasterblaster: libprasterblaster.so $(PRBOBJS) $(GCTPDIR)/libgctp_cpp.a
	$(CXX) $(LCFLAGS) $(LINCLUDES) src/driver.cpp -o prasterblaster $(LLIBS) -lprasterblaster

.PHONY:  clean

clean:
	rm -f libprasterblaster.a libprasterblaster.so prasterblaster src/*.o src/gctp_cpp/*.o src/gctp_cpp/*.a
