TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    celestialbody.cpp \
    solarsystem.cpp \
    integrator.cpp

HEADERS += \
    celestialbody.h \
    solarsystem.h \
    integrator.h

