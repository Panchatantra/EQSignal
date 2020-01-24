#-------------------------------------------------
#
# Project created by QtCreator 2015-07-10T12:16:11
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = EQSignal
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    pulse.cpp \
    qcustomplot.cpp \
    EQSignal.cpp \
    Spectra.cpp \
    eqtablewidget.cpp \
    genwavedialog.cpp \
    spectradefinewidget.cpp

HEADERS  += mainwindow.h \
    pulse.h \
    qcustomplot.h \
    EQSignal.h \
    eqs.h \
    Spectra.h \
    eqtablewidget.h \
    genwavedialog.h \
    spectradefinewidget.h

unix: QMAKE_CXXFLAGS += -std=c++11 -fopenmp

FORMS    += mainwindow.ui \
    pulse.ui

RESOURCES += mainwindow.qrc

TRANSLATIONS += trans.ts

unix:!macx: LIBS += -leqs -lbwf -lgomp -lpthread

macx|win32: LIBS += -L$$PWD/lib/ -leqs -lbwf

macx: ICON = icons/eqs.icns

win32: RC_FILE = EQSignal.rc

