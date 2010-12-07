SRCDIR = src/
PRBDIR = src/
GCTPDIR = src/gctp_cpp/
OBJDIR = src/objs

LCFLAGS = $(CFLAGS) -g -Wall -W -O2 -L$(GCTPDIR)
LINCLUDES = $(INCLUDES) -Isrc -I$(GCTPDIR)
LLIBS = $(LIBS) -lmpi -lgdal -lgctp_cpp

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
               driver.cpp \
               rasterinfo.cpp \
               reprojector.cpp \
               tinystr.cpp \
               tinyxml.cpp \
               tinyxmlerror.cpp \
               tinyxmlparser.cpp \
               rasterxml.cpp \
               resampler.cpp


GCTPSOURCES = $(patsubst %,$(GCTPDIR)/%,$(_GCTPSOURCES))
PRBSOURCES = $(patsubst %,$(PRBDIR)/%,$(_PRBSOURCES))
GCTPOBJS = $(patsubst %.cpp,%.o,$(GCTPSOURCES))
PRBOBJS = $(patsubst %.cpp,%.o,$(PRBSOURCES))

all: prasterblaster 

.cpp.o:
	$(CXX) $(LCFLAGS) $(LINCLUDES) -c $< -o $@

$(GCTPDIR)/libgctp_cpp.a: $(GCTPOBJS)
	ar rcs $(GCTPDIR)/libgctp_cpp.a $(GCTPOBJS)

prasterblaster: $(PRBOBJS) $(GCTPDIR)/libgctp_cpp.a
	$(CXX) $(LCFLAGS) $(LINCLUDES)  -o prasterblaster $(PRBOBJS) $(LLIBS)

.PHONY:  clean

clean:
	rm -f prasterblaster src/*.o src/gctp_cpp/*.o src/gctp_cpp/*.a
