#-------------------------------------------------
#
# Project created by QtCreator 2012-05-08T06:49:43
#
#-------------------------------------------------

QT       -= core gui

TARGET = librasterblaster
TEMPLATE = lib

DEFINES += LIBRASTERBLASTER_LIBRARY

SOURCES += projectedraster.cpp \
           rasterchunk.cpp \
           reprojector.cpp \
           resampler.cpp \

HEADERS += projectedraster.hh \
           rasterchunk.hh \
           reprojector.hh \
           resampler.hh \
           sharedptr.hh

QMAKE_CXXFLAGS += -DHAVE_TR1_SHARED_PTR

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
