######################################################################
# Automatically generated by qmake (2.01a) Tue Mar 13 16:34:47 2012
######################################################################


QMAKEFEATURES += /usr/lib64/qt4/features6
CONFIG += qwt
LIBS += -lqwt6


TEMPLATE = app
TARGET = solver
OBJECTS_DIR = ./tmp
DEPENDPATH += . \
              gui \
              material-data/can
INCLUDEPATH += . \
               /usr/include/qwt6 \
               material-data/can \
               gui

# Input
HEADERS += fd-solver.h \
           gui/solver.h \
           material-data/can/can.h \
           material-data/can/datafile.h
FORMS += gui/solver.ui
SOURCES += domain1d.c \
           node1d.c \
           output.c \
           update-functions.c \
           gui/main.cpp \
           gui/solver.cpp \
           material-data/can/can.c \
           material-data/can/datafile.c