#include "spectradefinewidget.h"

SpectraDefineWidget::SpectraDefineWidget(QWidget *parent) : QWidget(parent)
{
    wl = new QVBoxLayout(this);
    bl = new QHBoxLayout;

    this->setLayout(wl);

    readButton  = new QPushButton("Read",this);
    applyButton = new QPushButton("Apply",this);

    label_np = new QLabel("NP:",this);

    int NP = 30;
    spinBox_np = new QSpinBox(this);
    spinBox_np->setValue(NP);

    bl->addWidget(readButton);
    bl->addWidget(applyButton);
    bl->addWidget(label_np);
    bl->addWidget(spinBox_np);

    dataTable = new EQTableWidget(this);

    dataTable->setColumnCount(2);
    dataTable->setRowCount(NP);

    wl->addLayout(bl);
    wl->addWidget(dataTable);

    connect(spinBox_np,SIGNAL(valueChanged(int)),dataTable,SLOT(setRowNumber(int)));
    connect(dataTable,SIGNAL(rowNumberChanged(int)),spinBox_np,SLOT(setValue(int)));

}

PeriodsDefineWidget::PeriodsDefineWidget(QWidget *parent) :
    SpectraDefineWidget(parent)
{
    dataTable->setColumnCount(1);

    QStringList h;
    h << "Period";
    dataTable->setHorizontalHeaderLabels(h);
}

SPTDefineWidget::SPTDefineWidget(QWidget *parent) :
    SpectraDefineWidget(parent)
{
    dataTable->setColumnCount(2);

    QStringList h;
    h << "Period" << "SPT";
    dataTable->setHorizontalHeaderLabels(h);
}
