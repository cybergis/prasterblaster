#-------------------------------------------------
#
# Project created by QtCreator 2012-05-08T06:49:43
#
#-------------------------------------------------

QT       -= core gui

TARGET = librasterblaster
TEMPLATE = lib

DEFINES += LIBRASTERBLASTER_LIBRARY
DEPENDPATH += . \
              gctp_cpp \
              gctp_cpp/python \
              gtest/codegear \
              gtest/samples \
              gtest/src \
              gtest/test \
              gtest/fused-src/gtest \
              gtest/include/gtest \
              gtest/include/gtest/internal \
              gtest/xcode/Samples/FrameworkSample
INCLUDEPATH += . \
               gctp_cpp \
               gtest/samples \
               gtest/test \
               gtest/xcode/Samples/FrameworkSample

HEADERS += projectedraster.hh \
           quadtree.hh \
           rasterchunk.hh \
           reprojector.hh \
           resampler.hh \
           sharedptr.hh \
           gctp_cpp/alaskaconformal.h \
           gctp_cpp/albersConEqArea.h \
           gctp_cpp/azequidistant.h \
           gctp_cpp/constants.h \
           gctp_cpp/coordinate.h \
           gctp_cpp/cproj.h \
           gctp_cpp/equidistantc.h \
           gctp_cpp/equirectangular.h \
           gctp_cpp/genvertnsp.h \
           gctp_cpp/gnomonic.h \
           gctp_cpp/goodeh.h \
           gctp_cpp/hammer.h \
           gctp_cpp/hotinobmerc.h \
           gctp_cpp/intmollweide.h \
           gctp_cpp/lambertazimuth.h \
           gctp_cpp/lambertcc.h \
           gctp_cpp/mercator.h \
           gctp_cpp/miller.h \
           gctp_cpp/mollweide.h \
           gctp_cpp/oblatedeqarea.h \
           gctp_cpp/orthographic.h \
           gctp_cpp/polarstereo.h \
           gctp_cpp/polyconic.h \
           gctp_cpp/projection.h \
           gctp_cpp/robinson.h \
           gctp_cpp/sinusoidal.h \
           gctp_cpp/spaceobmerc.h \
           gctp_cpp/stereo.h \
           gctp_cpp/transformer.h \
           gctp_cpp/transversemercator.h \
           gctp_cpp/util.h \
           gctp_cpp/utm.h \
           gctp_cpp/vandergrinten.h \
           gctp_cpp/wagneriv.h \
           gctp_cpp/wagnervii.h 

SOURCES +=            projectedraster.cpp \
           quadtree.cpp \
           rasterchunk.cpp \
           reprojector.cpp \
           resampler.cpp \
           gctp_cpp/alaskaconformal.cpp \
           gctp_cpp/albersConEqArea.cpp \
           gctp_cpp/azequidistant.cpp \
           gctp_cpp/constants.cpp \
           gctp_cpp/coordinate.cpp \
           gctp_cpp/cproj.cpp \
           gctp_cpp/equidistantc.cpp \
           gctp_cpp/equirectangular.cpp \
           gctp_cpp/genvertnsp.cpp \
           gctp_cpp/gnomonic.cpp \
           gctp_cpp/goodeh.cpp \
           gctp_cpp/hammer.cpp \
           gctp_cpp/hotinobmerc.cpp \
           gctp_cpp/intmollweide.cpp \
           gctp_cpp/lambertazimuth.cpp \
           gctp_cpp/lambertcc.cpp \
           gctp_cpp/main.cpp \
           gctp_cpp/mercator.cpp \
           gctp_cpp/miller.cpp \
           gctp_cpp/Mollweide.cpp \
           gctp_cpp/oblatedeqarea.cpp \
           gctp_cpp/orthographic.cpp \
           gctp_cpp/polarstereo.cpp \
           gctp_cpp/polyconic.cpp \
           gctp_cpp/projection.cpp \
           gctp_cpp/robinson.cpp \
           gctp_cpp/sinusoidal.cpp \
           gctp_cpp/spaceobmerc.cpp \
           gctp_cpp/stereo.cpp \
           gctp_cpp/transformer.cpp \
           gctp_cpp/transversemercator.cpp \
           gctp_cpp/util.cpp \
           gctp_cpp/utm.cpp \
           gctp_cpp/vandergrinten.cpp \
           gctp_cpp/wagneriv.cpp \
           gctp_cpp/wagnervii.cpp \

QMAKE_CXXFLAGS += -DHAVE_TR1_SHARED_PTR
HAVE_CONFIG_H = 0

symbian {
    MMP_RULES += EXPORTUNFROZEN
    TARGET.UID3 = 0xE0A22788
    TARGET.CAPABILITY = 
    TARGET.EPOCALLOWDLLDATA = 1
    addFiles.sources = librasterblaster.dll
    addFiles.path = !:/sys/bin
    DEPLOYMENT += addFiles
}

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
