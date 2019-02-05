TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -larmadillo

SOURCES += main.cpp \
    LU_bench.cpp

