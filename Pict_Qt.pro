QT += core
QT -= gui
QT += widgets

CONFIG += c++11

TARGET = Pict_Qt

CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    haar.cpp \
    subimage.cpp \
    src/cluster.c \
    mainwindow.cpp \
    clustsigma.cpp

INCLUDEPATH += c:/boost/boost_1_62_0
LIBS += "-LC:/boost_lib/"

HEADERS += \
    haar.hpp \
    subimage.h \
    src/cluster.h \
    mainwindow.h \
    clustsigma.h

FORMS += \
    mainwindow.ui

RESOURCES += \
    resources.qrc
