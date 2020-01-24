#include "spectradefinewidget.h"

SpectraDefineWidget::SpectraDefineWidget(QWidget *parent) : QWidget(parent)
{
    wl = new QVBoxLayout(this);
    bl = new QHBoxLayout;
    drl = new QVBoxLayout;

    this->setLayout(wl);

    readButton  = new QPushButton(tr("From File"),this);
    applyButton = new QPushButton(tr("Apply"),this);

	label_np = new QLabel(tr("NP"), this);
	label_dr = new QLabel(tr("Damping Reduction Relationship"), this);
	label_table = new QLabel(tr("Define 5%-Damping Spectrum"), this);

    int NP = 60;
    spinBox_np = new QSpinBox(this);
    spinBox_np->setValue(NP);
	spinBox_np->setMaximum(10000);

	comBox_dr = new QComboBox(this);
	comBox_dr->addItem(tr("GB 50011-2010"));
	comBox_dr->addItem(tr("FEMA 440"));

    bl->addWidget(readButton);
    bl->addWidget(applyButton);
    bl->addWidget(label_np);
    bl->addWidget(spinBox_np);

	drl->addWidget(label_dr);
	drl->addWidget(comBox_dr);
	drl->addWidget(label_table);

    dataTable = new EQTableWidget(this);

    dataTable->setColumnCount(2);
    dataTable->setRowCount(NP);

	wl->addLayout(bl);
	wl->addLayout(drl);
    wl->addWidget(dataTable);

    connect(spinBox_np,SIGNAL(valueChanged(int)),dataTable,SLOT(setRowNumber(int)));
	connect(dataTable, SIGNAL(rowNumberChanged(int)), spinBox_np, SLOT(setValue(int)));

}

PeriodsDefineWidget::PeriodsDefineWidget(QWidget *parent) :
    SpectraDefineWidget(parent)
{
    dataTable->setColumnCount(1);

    QStringList h;
    h << tr("Period");
    dataTable->setHorizontalHeaderLabels(h);

}

void PeriodsDefineWidget::readFromFile()
{
	QString fileName, c;
	QStringList filters, clist;
	QFileDialog fdialog(this);

	//filters << "Supported file (*.txt *.TXT *.csv *.CSV)"
	//	<< "Text file (*.txt *.TXT)"
	//	<< "CSV file (*.csv *.CSV)";

	filters << "Text file (*.txt *.TXT)";

	fdialog.setNameFilters(filters);
	fdialog.setWindowTitle(tr("From File"));

	if (fdialog.exec()) fileName = fdialog.selectedFiles()[0];

	if (fileName.isEmpty()) return;

	QFile file(fileName);
	file.open(QIODevice::ReadOnly);
	QTextStream ts(&file);

	c = ts.readAll();
    clist = c.split(QRegExp("\\s+"));

	int np = clist.count();
	setNP(np);
	setPeriods(clist);

}

SPTDefineWidget::SPTDefineWidget(QWidget *parent) :
    SpectraDefineWidget(parent)
{
    dataTable->setColumnCount(2);

    QStringList h;
    h << tr("Period") << tr("SPT");
    dataTable->setHorizontalHeaderLabels(h);
}

void SPTDefineWidget::readFromFile()
{
	QString fileName, c;
	QStringList filters, clist, plist, sptlist;
	QFileDialog fdialog(this);

	//filters << "Supported file (*.txt *.TXT *.csv *.CSV)"
	//	<< "Text file (*.txt *.TXT)"
	//	<< "CSV file (*.csv *.CSV)";

	filters  << "Text file (*.txt *.TXT)";

	fdialog.setNameFilters(filters);
	fdialog.setWindowTitle(tr("From File"));

	if (fdialog.exec()) fileName = fdialog.selectedFiles()[0];

	if (fileName.isEmpty()) return;

	QFile file(fileName);
	file.open(QIODevice::ReadOnly);
	QTextStream ts(&file);

	c = ts.readAll();
	clist = c.split(QRegExp("\\s+"));

	int np = clist.count()/2;
	setNP(np);

	for (int i = 0; i < np; i++)
	{
		plist << clist[2 * i];
		sptlist << clist[2 * i + 1];
	}

	setPeriods(plist);
	setSPT(sptlist);
}
