TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

SOURCES += main.cpp \
    montecarlostep.cpp \
    meanvalue.cpp \
    particle.cpp \
    system.cpp \
    numericaldiff.cpp

HEADERS += \
    montecarlostep.h \
    meanvalue.h \
    particle.h \
    system.h \
    numericaldiff.h

