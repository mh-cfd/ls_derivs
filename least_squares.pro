TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    main.cpp \
    leastsquaressolver.cpp

HEADERS += \
    leastsquaressolver.h

QMAKE_CXXFLAGS_RELEASE += -O3


QMAKE_LFLAGS += -O3


#LIBS += -lOpenGL32 -lGLU32  -lm
#LIBS += -L$$PWD/my_lib -lglut32


LIBS += -lGL -lGLU -lglut  -lm

