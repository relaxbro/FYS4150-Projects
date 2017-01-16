TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/gaussiandeviate.cpp\
    src/particle.cpp \
    src/barneshut.cpp \
    src/octree.cpp \
    src/clustersystem.cpp

HEADERS += \
    src/particle.h \
    src/barneshut.h \
    src/octree.h \
    src/clustersystem.h

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS_RELEASE += -fopenmp
    QMAKE_LFLAGS_RELEASE += -fopenmp
}
