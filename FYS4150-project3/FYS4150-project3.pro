TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    src/solarsystem.cpp \
    src/solarsystemobject.cpp \
    src/verlet.cpp \
    src/rungekutta4.cpp

HEADERS += \
    src/solarsystem.h \
    src/solarsystemobject.h \
    src/verlet.h \
    src/rungekutta4.h

#release {
#    DEFINES += ARMA_NO_DEBUG
#    QMAKE_CXXFLAGS_RELEASE -= O2
#    QMAKE_CXXFLAGS_RELEASE += O3
#}

#OTHER_FILES +=

copydata.commands = $(COPY_DIR) $$PWD/initialdata $$OUT_PWD
first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata

OTHER_FILES += \
    initialdata/bodies.txt
