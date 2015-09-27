#ifndef SPECTRADEFINEWIDGET
#define SPECTRADEFINEWIDGET

#include <QtCore>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSpinBox>
#include <QVBoxLayout>
#include <QWidget>
#include "eqtablewidget.h"

class SpectraDefineWidget : public QWidget
{
    Q_OBJECT
public:
    SpectraDefineWidget(QWidget *parent=0);

    QVBoxLayout *wl;
    QHBoxLayout *bl;

    QPushButton *readButton;
    QPushButton *applyButton;

    QLabel *label_np;
    QSpinBox *spinBox_np;

    EQTableWidget *dataTable;

};

class PeriodsDefineWidget : public SpectraDefineWidget
{
    Q_OBJECT
public:
    PeriodsDefineWidget(QWidget *parent=0);
    double *getPeriods() {return dataTable->getColumnData(0);}

};

class SPTDefineWidget : public SpectraDefineWidget
{
    Q_OBJECT
public:
    SPTDefineWidget(QWidget *parent=0);

    double *getPeriods() {return dataTable->getColumnData(0);}
    double *getSPT() {return dataTable->getColumnData(1);}

};

#endif // SPECTRADEFINEWIDGET

