TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_LIBS += -ltiff

SOURCES += main.cpp \
    enhancement.cpp \
    bmp.c \
    image.cpp

HEADERS += \
    enhancement.hpp
