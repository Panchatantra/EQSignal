#ifndef SPECTRADEFINEWIDGET
#define SPECTRADEFINEWIDGET

#include <QtCore>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSpinBox>
#include <QComboBox>
#include <QVBoxLayout>
#include <QWidget>
#include <QFileDialog>
#include "eqtablewidget.h"

class SpectraDefineWidget : public QWidget
{
    Q_OBJECT

public:

    SpectraDefineWidget(QWidget *parent=0);

    QVBoxLayout *wl;
    QHBoxLayout *bl;
    QVBoxLayout *drl;

    QPushButton *readButton;
    QPushButton *applyButton;

	QLabel *label_np;
	QLabel *label_dr;
	QLabel *label_table;
    QSpinBox *spinBox_np;
    QComboBox *comBox_dr;

    EQTableWidget *dataTable;

    void setNP(int np) {spinBox_np->setValue(np);}
	void setColumn(int col, double *data) { dataTable->setColumn(col, data); }
	void setColumn(int col, QStringList data) { dataTable->setColumn(col, data); }
	void setPeriods(double *data) { setColumn(0, data); }
	void setPeriods(QStringList data) { setColumn(0, data); }

	virtual void readFromFile() {};

};

class PeriodsDefineWidget : public SpectraDefineWidget
{
    Q_OBJECT

public:

    PeriodsDefineWidget(QWidget *parent=0);
    double *getPeriods() {return dataTable->getColumnData(0);}

	void readFromFile();

};

class SPTDefineWidget : public SpectraDefineWidget
{
    Q_OBJECT

public:

    SPTDefineWidget(QWidget *parent=0);

    double *getPeriods() {return dataTable->getColumnData(0);}
    double *getSPT() {return dataTable->getColumnData(1);}
	void setSPT(double *data) { setColumn(1, data); }
	void setSPT(QStringList data) { setColumn(1, data); }

	void readFromFile();

};

#endif // SPECTRADEFINEWIDGET

