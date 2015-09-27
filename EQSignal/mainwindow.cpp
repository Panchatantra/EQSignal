#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "qcustomplot.h"
#include "eqs.h"
#include <time.h>
#include <algorithm>
#include <fstream>
#include <QVector>
#include <QMessageBox>
#include <QProgressDialog>
#include <QJsonDocument>
#include <QJsonObject>

using namespace std;

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent),
    ui(new Ui::MainWindow),
    gwd(new GenWaveDialog),
    eqs(new EQSignal),
    eqs0(new EQSignal),
    pdw(new PeriodsDefineWidget),
    sptw(new SPTDefineWidget),
    wth(new QWidget),
    wres(new QWidget),
    wsp(new QWidget)
{
    ui->setupUi(this);

    initViewTH();
    initViewSPA();
    initViewSP();
    initViewFFT();
    initViewRES();
    initViewEnergy();
    initViewHyst();
    initTable();
    setupConnections();

    workDir = QDir::current();
    saveAccWithTime = true;
    normOnRead = false;
    FIT_SPA = false;

    this->readConfig();
    ui->Paras->setCurrentIndex(0);
    ui->tabWidget->setCurrentIndex(0);

    ui->progressBar->hide();

    #ifndef Q_OS_MAC
    ui->dockToolBox->setFloating(true);
    ui->dockXScale->setFloating(true);
    ui->dockToolBox->move(10, 0);
    ui->dockToolBox->resize(300,640);
    ui->dockXScale->move(10, 660);
    ui->dockXScale->resize(300,90);
    #endif

    #ifdef Q_OS_MAC
    ui->menuBar->setParent(0);
	#endif

}

MainWindow::~MainWindow()
{
	delete ui;
	delete eqs;
    delete eqs0;
}

void MainWindow::setupConnections()
{
    connect(ui->actionGenWave, &QAction::triggered,
            gwd, &GenWaveDialog::show);
    connect(gwd, &GenWaveDialog::accepted,
            this, &MainWindow::genWave);

    connect(ui->Open, &QToolButton::clicked,
            this, &MainWindow::on_actionOpen_triggered);
    connect(ui->Confirm, &QPushButton::clicked,
            this, &MainWindow::confirm);
    connect(ui->Confirm_Pre, &QPushButton::clicked,
            this, &MainWindow::confirm);

    connect(ui->XS_LIN_SP, &QRadioButton::clicked,
            this, &MainWindow::SPXScaleChanged);
    connect(ui->XS_LOG_SP, &QRadioButton::clicked,
            this, &MainWindow::SPXScaleChanged);
    connect(ui->Amp, &QRadioButton::clicked,
            this, &MainWindow::plotFFT);
    connect(ui->Ang, &QRadioButton::clicked,
            this, &MainWindow::plotFFT);
    connect(ui->DAng, &QRadioButton::clicked,
            this, &MainWindow::plotFFT);
    connect(ui->IAON, &QCheckBox::clicked,
            this, &MainWindow::plotTH);
    connect(ui->INPUTON, &QCheckBox::clicked,
            this, &MainWindow::plotRES);
    connect(dr,static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,&MainWindow::fillSPTable);
    connect(pdw->applyButton, &QPushButton::clicked,
            this, &MainWindow::pdw_applyButton_clicked);
    connect(sptw->applyButton, &QPushButton::clicked,
            this, &MainWindow::sptw_applyButton_clicked);
}

void MainWindow::readtxt(const char *filename, double DT)
{
	eqs->readtxt(filename, DT);
	eqs0->readtxt(filename, DT);

	int N = eqs->getN();
	ui->Ind1->setValue(0);
	ui->Ind2->setValue(N - 1);
	ui->Ind2->setMaximum(N - 1);
	ui->TS->setValue(eqs->getDt()*2.5);

    plotTH();
}

void MainWindow::readtxt(QString filename, double DT)
{
    eqs->readtxt(filename, DT, normOnRead);
    eqs0->readtxt(filename, DT, normOnRead);

	int N = eqs->getN();
	ui->Ind1->setValue(0);
	ui->Ind2->setValue(N - 1);
	ui->Ind2->setMaximum(N - 1);
	ui->TS->setValue(eqs->getDt()*2.5);

    plotTH();
}

void MainWindow::readnga(const char *filename)
{
    eqs->readnga(filename, normOnRead);
    eqs0->readnga(filename, normOnRead);

    int N = eqs->getN();
    ui->Ind1->setValue(0);
    ui->Ind2->setValue(N - 1);
    ui->Ind2->setMaximum(N - 1);
	ui->TS->setValue(eqs->getDt()*2.5);

    this->plotTH();
}

void MainWindow::readnga(QString filename)
{
    eqs->readnga(filename, normOnRead);
    eqs0->readnga(filename, normOnRead);

    int N = eqs->getN();
    ui->Ind1->setValue(0);
    ui->Ind2->setValue(N - 1);
    ui->Ind2->setMaximum(N - 1);
	ui->TS->setValue(eqs->getDt()*2.5);

    this->plotTH();
}

QPen MainWindow::autoPen(int i)
{
    int nc = 6;
    if (i >= nc) i = i % nc;

    QPen pen;

    switch (i)
    {
    case 0:
        pen = QPen(QColor(  0,  0,255)); // blue
        break;
    case 1:
        pen = QPen(QColor(255,  0,  0)); // red
        break;
    case 2:
        pen = QPen(QColor(  0,255,  0)); // green
        break;
    case 3:
        pen = QPen(QColor(255,  0,255)); // Magenta
        break;
    case 4:
        pen = QPen(QColor(  0,255,255)); // cyan
        break;
    case 5:
        pen = QPen(QColor(255,255,  0)); // yellow
        break;
    default:
        break;
    }

    pen.setWidthF(1.5);

    return pen;
}

QColor MainWindow::reverseColor(QColor c)
{
    int r = 255 - c.red();
    int g = 255 - c.green();
    int b = 255 - c.blue();

    return QColor(r,g,b);
}

void MainWindow::initViewSP()
{
    QCustomPlot *qplot = ui->ViewSP;
    qplot->plotLayout()->clear();

	QCPAxisRect *plotSPA = new QCPAxisRect(qplot);
	QCPAxisRect *plotSPV = new QCPAxisRect(qplot);
    QCPAxisRect *plotSPD = new QCPAxisRect(qplot);
	QCPAxisRect *plotSPE = new QCPAxisRect(qplot);

	qplot->plotLayout()->addElement(0, 0, plotSPA);
    qplot->plotLayout()->addElement(0, 1, plotSPV);
	qplot->plotLayout()->addElement(1, 0, plotSPD);
	qplot->plotLayout()->addElement(1, 1, plotSPE);

	QCPMarginGroup *marginGroup = new QCPMarginGroup(qplot);

    plotSPA->setMarginGroup(QCP::msAll, marginGroup);
    plotSPV->setMarginGroup(QCP::msAll, marginGroup);
    plotSPD->setMarginGroup(QCP::msAll, marginGroup);
    plotSPE->setMarginGroup(QCP::msAll, marginGroup);

	foreach(QCPAxisRect *rect, qplot->axisRects())
	{
        rect->axis(QCPAxis::atBottom)->setRange(0.01, 10.0);
        if (ui->XS_LOG_SP->isChecked())
        {
            rect->axis(QCPAxis::atBottom)->setScaleType(QCPAxis::stLogarithmic);
        }
        else if (ui->XS_LIN_SP->isChecked())
        {
            rect->axis(QCPAxis::atBottom)->setScaleType(QCPAxis::stLinear);
        }
        rect->axis(QCPAxis::atBottom)->grid()->setSubGridVisible(true);
		rect->axis(QCPAxis::atBottom)->setLabel(tr("Period"));
        rect->axis(QCPAxis::atLeft)->setRange(0.0, 2.5);

        foreach(QCPAxis *axis, rect->axes())
		{
			axis->setLayer("axes");
			axis->grid()->setLayer("grid");
		}
	}

	plotSPA->axis(QCPAxis::atLeft)->setLabel(tr("SPA"));
	plotSPV->axis(QCPAxis::atLeft)->setLabel(tr("SPV"));
	plotSPD->axis(QCPAxis::atLeft)->setLabel(tr("SPD"));
	plotSPE->axis(QCPAxis::atLeft)->setLabel(tr("SPE"));

}

void MainWindow::initViewFFT()
{
    QCustomPlot *qplot = ui->ViewFFT;

    qplot->setInteractions(QCP::iRangeZoom|QCP::iRangeDrag);

    qplot->axisRect()->setRangeZoom(Qt::Horizontal);
    qplot->axisRect()->setRangeDrag(Qt::Horizontal);

    qplot->xAxis->setLabel(tr("Frequency (Hz)"));
    qplot->yAxis->setLabel(tr("Fourier Amplitude"));
    qplot->xAxis->setRange(0.0,50.0);
    qplot->yAxis->setRange(0.0,1.0);

    qplot->xAxis->grid()->setSubGridVisible(true);

}

void MainWindow::initTable()
{

    dataTableTH = new EQTableWidget(wth);
    QVBoxLayout *lth = new QVBoxLayout(wth);

    QStringList HHeader;
    HHeader << tr("Time") << tr("Acceleration") << tr("Velocity") << tr("Displacement");

    wth->resize(480,600);
    wth->setWindowIcon(QIcon(":/MainWindow/icons/ViewData.png"));
    wth->setWindowTitle(tr("Time History Data"));

    dataTableTH->setColumnCount(4);
    dataTableTH->setHorizontalHeaderLabels(HHeader);

    lth->addWidget(dataTableTH);

    dataTableRES = new EQTableWidget(wres);
    QVBoxLayout *lres = new QVBoxLayout(wres);

    wres->resize(480,600);
    wres->setWindowIcon(QIcon(":/MainWindow/icons/ViewData.png"));
    wres->setWindowTitle(tr("Time History Data"));

    dataTableTH->setColumnCount(4);
    dataTableTH->setHorizontalHeaderLabels(HHeader);

    lres->addWidget(dataTableRES);

    HHeader.clear();
    HHeader << tr("Period") << tr("SPA") << tr("SPV") << tr("SPD") << tr("SPE") << tr("SPT");

    QHBoxLayout *ldr = new QHBoxLayout;
    QVBoxLayout *lsp = new QVBoxLayout(wsp);
    QLabel *l = new QLabel(tr("Damping Ratio"),wsp);
    dr = new QComboBox(wsp);

    dataTableSP = new EQTableWidget(wsp);

    wsp->resize(640,600);
    wsp->setWindowIcon(QIcon(":/MainWindow/icons/ViewData.png"));
    wsp->setWindowTitle(tr("Response Spectra Data"));
    dataTableSP->setColumnCount(6);
    dataTableSP->setHorizontalHeaderLabels(HHeader);

    dr->addItem("0.05");
    dr->addItem("0.02");

    l->setAlignment(Qt::AlignRight);
    ldr->addWidget(l);
    ldr->addWidget(dr);
    lsp->addLayout(ldr);
    lsp->addWidget(dataTableSP);

    pdw->setWindowIcon(QIcon(":/MainWindow/icons/eqs.ico"));
    sptw->setWindowIcon(QIcon(":/MainWindow/icons/eqs.ico"));

    //    wsp->show();
}

void MainWindow::clearView(QCustomPlot *qplot)
{
    qplot->clearGraphs();
}

void MainWindow::readConfig()
{
    QString f = QDir::home().filePath(".eqsignal");
    if (!QFile::exists(f)) return;

    QFile file(f);
    file.open(QIODevice::ReadOnly);
    QTextStream ts(&file);

    QString conf = ts.readAll();
    QJsonObject json = QJsonDocument::fromJson(conf.toUtf8()).object();

    workDir = QDir(json["work_dir"].toString());
    lastFile = json["last_file"].toString();

    saveAccWithTime = json["save_acc_with_time"].toBool();
    normOnRead = json["norm_on_read"].toBool();

    ui->URL->setText(lastFile);

    file.close();

}

void MainWindow::writeConfig()
{
    QString f = QDir::home().filePath(".eqsignal");
    if (!QFile::exists(f)) return;

    QFile file(f);
    file.open(QIODevice::WriteOnly);

    QJsonObject json;
    json.insert("work_dir",QJsonValue(workDir.path()));
    json.insert("save_acc_with_time",QJsonValue(saveAccWithTime));
    json.insert("norm_on_read",QJsonValue(normOnRead));
    json.insert("last_file",QJsonValue(lastFile));

    QJsonDocument jdoc(json);

    file.write(jdoc.toJson());
    file.close();

}

void MainWindow::initViewTH()
{
    QCustomPlot *qplot = ui->ViewTH;
    qplot->plotLayout()->clear();
    QCPAxisRect *plotAcc = new QCPAxisRect(qplot);
    QCPAxisRect *plotVel = new QCPAxisRect(qplot);
    QCPAxisRect *plotDsp = new QCPAxisRect(qplot);

    qplot->plotLayout()->addElement(0, 0, plotAcc);
    qplot->plotLayout()->addElement(1, 0, plotVel);
    qplot->plotLayout()->addElement(2, 0, plotDsp);

    QCPMarginGroup *marginGroup = new QCPMarginGroup(qplot);

    plotAcc->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);
    plotVel->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);
    plotDsp->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);

    foreach(QCPAxisRect *rect, qplot->axisRects())
    {
        foreach(QCPAxis *axis, rect->axes())
        {
            axis->setLayer("axes");
            axis->grid()->setLayer("grid");
        }
    }

    plotAcc->axis(QCPAxis::atBottom)->setRange(0.0, 40.0);
    plotVel->axis(QCPAxis::atBottom)->setRange(0.0, 40.0);
    plotDsp->axis(QCPAxis::atBottom)->setRange(0.0, 40.0);

    plotAcc->axis(QCPAxis::atLeft)->setRange(-1.0, 1.0);
    plotVel->axis(QCPAxis::atLeft)->setRange(-1.0, 1.0);
    plotDsp->axis(QCPAxis::atLeft)->setRange(-1.0, 1.0);

    plotAcc->axis(QCPAxis::atLeft)->setLabel(tr("Acceleration"));
    plotVel->axis(QCPAxis::atLeft)->setLabel(tr("Velocity"));
    plotDsp->axis(QCPAxis::atLeft)->setLabel(tr("Displacement"));
    plotDsp->axis(QCPAxis::atBottom)->setLabel(tr("Time"));

//    qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atLeft));
//    qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atLeft));
//    qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atLeft));

//    qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atLeft));
//    qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atLeft));
//    qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atLeft));

//    qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atLeft));
//    qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atLeft));
//    qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atLeft));

//    QPen penRaw = QPen(Qt::gray);
//    qplot->graph(0)->setPen(penRaw);
//    qplot->graph(1)->setPen(penRaw);
//    qplot->graph(2)->setPen(penRaw);

//    QPen penTar = QPen(QColor(46, 139, 87));
//    qplot->graph(6)->setPen(penTar);
//    qplot->graph(7)->setPen(penTar);
    //    qplot->graph(8)->setPen(penTar);
}

void MainWindow::initViewRES()
{
    QCustomPlot *qplot = ui->ViewRES;
    qplot->plotLayout()->clear();
    QCPAxisRect *plotAcc = new QCPAxisRect(qplot);
    QCPAxisRect *plotVel = new QCPAxisRect(qplot);
    QCPAxisRect *plotDsp = new QCPAxisRect(qplot);

    qplot->plotLayout()->addElement(0, 0, plotAcc);
    qplot->plotLayout()->addElement(1, 0, plotVel);
    qplot->plotLayout()->addElement(2, 0, plotDsp);

    QCPMarginGroup *marginGroup = new QCPMarginGroup(qplot);

    plotAcc->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);
    plotVel->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);
    plotDsp->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);

    foreach(QCPAxisRect *rect, qplot->axisRects())
    {
        foreach(QCPAxis *axis, rect->axes())
        {
            axis->setLayer("axes");
            axis->grid()->setLayer("grid");
        }
    }

    plotAcc->axis(QCPAxis::atBottom)->setRange(0.0, 40.0);
    plotVel->axis(QCPAxis::atBottom)->setRange(0.0, 40.0);
    plotDsp->axis(QCPAxis::atBottom)->setRange(0.0, 40.0);

    plotAcc->axis(QCPAxis::atLeft)->setRange(-1.0, 1.0);
    plotVel->axis(QCPAxis::atLeft)->setRange(-1.0, 1.0);
    plotDsp->axis(QCPAxis::atLeft)->setRange(-1.0, 1.0);

    plotAcc->axis(QCPAxis::atLeft)->setLabel(tr("Acceleration"));
    plotVel->axis(QCPAxis::atLeft)->setLabel(tr("Velocity"));
    plotDsp->axis(QCPAxis::atLeft)->setLabel(tr("Displacement"));
    plotDsp->axis(QCPAxis::atBottom)->setLabel(tr("Time"));
}

void MainWindow::initViewHyst()
{
    QCustomPlot *qplot = ui->ViewHyst;

    qplot->legend->setVisible(false);
    qplot->legend->setBrush(QBrush(QColor(255,255,255,150)));

//    qplot->setInteractions(QCP::iRangeZoom|QCP::iRangeDrag);

//    qplot->axisRect()->setRangeZoom(Qt::Horizontal);
//    qplot->axisRect()->setRangeDrag(Qt::Horizontal);

    qplot->xAxis->setRange(-1.0, 1.0);
    qplot->yAxis->setRange(-1.0, 1.0);

    qplot->xAxis->setLabel(tr("Displacement"));
    qplot->yAxis->setLabel(tr("Force"));
}

void MainWindow::initViewSPA()
{
	QCustomPlot *qplot = ui->ViewSPA;

    qplot->legend->setVisible(false);
    qplot->legend->setBrush(QBrush(QColor(255,255,255,150)));

    qplot->setInteractions(QCP::iRangeZoom|QCP::iRangeDrag);

    qplot->axisRect()->setRangeZoom(Qt::Horizontal);
    qplot->axisRect()->setRangeDrag(Qt::Horizontal);

    if (ui->XS_LOG_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    }
    else if (ui->XS_LIN_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLinear);
    }

	qplot->xAxis->grid()->setSubGridVisible(true);
	qplot->xAxis->setRange(0.01, 10.0);
    qplot->yAxis->setRange(0.0, 3.0);

	qplot->xAxis->setLabel(tr("Period"));
    qplot->yAxis->setLabel(tr("Response Acceleration"));

}

void MainWindow::initViewEnergy()
{
    QCustomPlot *qplot = ui->ViewEnergy;

    qplot->legend->setVisible(false);
    qplot->legend->setBrush(QBrush(QColor(255,255,255,150)));

//    qplot->setInteractions(QCP::iRangeZoom|QCP::iRangeDrag);

//    qplot->axisRect()->setRangeZoom(Qt::Horizontal);
//    qplot->axisRect()->setRangeDrag(Qt::Horizontal);

    qplot->xAxis->setRange(0.0, 60.0);
    qplot->yAxis->setRange(0.0, 1.0);

    qplot->xAxis->setLabel(tr("Time"));
    qplot->yAxis->setLabel(tr("Energy"));

}

void MainWindow::saveSP()
{
	QString filename, sf;
	QStringList filters;
	QFileDialog fdialog(this);

	filters << "Text file (*.txt)" << "CSV file (*.csv)";

	fdialog.setAcceptMode(QFileDialog::AcceptSave);
	fdialog.setNameFilters(filters);
	fdialog.setWindowTitle(tr("Save"));
	fdialog.setDirectory(workDir);
	fdialog.selectFile(eqsName + "-SP");

	if (fdialog.exec()) {
		filename = fdialog.selectedFiles()[0];
		sf = fdialog.selectedNameFilter();
		if (sf == filters[0]) {
			eqs->savetxtsp(filename);
		}
		else if (sf == filters[1]) {
			eqs->savecsvsp(filename);
		}
	}
	
	//QString filename,autoname;
 //   autoname.append(getenv("EQS_WORK_DIR"));
 //   autoname.append(eqsName);
 //   autoname.append("-SP.txt");
 //   filename = qfile->getSaveFileName(this, "Save", autoname, "*.txt *.csv");

    //if (filename.endsWith(".txt")) eqs->savetxtsp(filename);
    //else if (filename.endsWith(".csv")) eqs->savecsvsp(filename);

}

void MainWindow::plotTH(bool changeTab)
{
    QCustomPlot *qplot = ui->ViewTH;
    qplot->clearGraphs();
    qplot->plotLayout()->clear();

    QCPAxisRect *plotAcc = new QCPAxisRect(qplot);
    QCPAxisRect *plotVel = new QCPAxisRect(qplot);
    QCPAxisRect *plotDsp = new QCPAxisRect(qplot);

    qplot->plotLayout()->addElement(0, 0, plotAcc);
    qplot->plotLayout()->addElement(1, 0, plotVel);
    qplot->plotLayout()->addElement(2, 0, plotDsp);

    QCPMarginGroup *marginGroup = new QCPMarginGroup(qplot);

    plotAcc->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);
    plotVel->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);
    plotDsp->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);

    foreach(QCPAxisRect *rect, qplot->axisRects())
    {
        rect->axis(QCPAxis::atBottom)->setRange(0.0, 40.0);
        rect->axis(QCPAxis::atLeft)->setRange(-1.0, 1.0);
        foreach(QCPAxis *axis, rect->axes())
        {
            axis->setLayer("axes");
            axis->grid()->setLayer("grid");
        }
    }

    plotAcc->axis(QCPAxis::atLeft)->setLabel(tr("Acceleration"));
    plotVel->axis(QCPAxis::atLeft)->setLabel(tr("Velocity"));
    plotDsp->axis(QCPAxis::atLeft)->setLabel(tr("Displacement"));
    plotDsp->axis(QCPAxis::atBottom)->setLabel(tr("Time"));

    qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atLeft));
    qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atLeft));
    qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atLeft));

    qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atLeft));
    qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atLeft));
    qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atLeft));

    qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atLeft));
    qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atLeft));
    qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atLeft));

    QPen penRaw = QPen(Qt::gray);
    qplot->graph(0)->setPen(penRaw);
    qplot->graph(1)->setPen(penRaw);
    qplot->graph(2)->setPen(penRaw);

    QPen penTar = QPen(QColor(46, 139, 87));
    qplot->graph(6)->setPen(penTar);
    qplot->graph(7)->setPen(penTar);
    qplot->graph(8)->setPen(penTar);

    QVector<double> Time = eqs->qGetT();
    bool e = false;

    if (ui->RawON->isChecked() || ui->RawONAlign->isChecked())
    {
        QVector<double> Acc0 = eqs0->qGetAcc();
        QVector<double> Vel0 = eqs0->qGetVel();
        QVector<double> Dsp0 = eqs0->qGetDisp();

        qplot->graph(0)->setData(Time, Acc0);
        qplot->graph(0)->rescaleAxes(e);
        qplot->graph(1)->setData(Time, Vel0);
        qplot->graph(1)->rescaleAxes(e);
        qplot->graph(2)->setData(Time, Dsp0);
        qplot->graph(2)->rescaleAxes(e);

        e = true;
    }

	QVector<double> Acc = eqs->qGetAcc();
	QVector<double> Vel = eqs->qGetVel();
	QVector<double> Dsp = eqs->qGetDisp();

	qplot->graph(3)->setData(Time, Acc);
	qplot->graph(3)->rescaleAxes(e);
	qplot->graph(4)->setData(Time, Vel);
	qplot->graph(4)->rescaleAxes(e);
	qplot->graph(5)->setData(Time, Dsp);
	qplot->graph(5)->rescaleAxes(e);

    if (ui->IAON->isChecked()) {
        eqs->calcAriasIntensity();
        QCPGraph *GIA = qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atRight));
        plotAcc->axis(QCPAxis::atRight)->setVisible(true);
        plotAcc->axis(QCPAxis::atRight)->setLabel(tr("Arias Intensity"));
        GIA->setData(Time,eqs->qGetIa());
        GIA->setPen(QPen(Qt::gray));
        GIA->setBrush(QBrush(QColor(0, 0, 255, 20)));
        GIA->rescaleValueAxis();

        QCPGraph *GIV = qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atRight));
        plotVel->axis(QCPAxis::atRight)->setVisible(true);
        plotVel->axis(QCPAxis::atRight)->setLabel(tr("Arias Intensity"));
        GIV->setData(Time,eqs->qGetIv());
        GIV->setPen(QPen(Qt::gray));
        GIV->setBrush(QBrush(QColor(0, 0, 255, 20)));
        GIV->rescaleValueAxis();

        QCPGraph *GID = qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atRight));
        plotDsp->axis(QCPAxis::atRight)->setVisible(true);
        plotDsp->axis(QCPAxis::atRight)->setLabel(tr("Arias Intensity"));
        GID->setData(Time,eqs->qGetId());
        GID->setPen(QPen(Qt::gray));
        GID->setBrush(QBrush(QColor(0, 0, 255, 20)));
        GID->rescaleValueAxis();

    }

	qplot->replot();
    if (changeTab) {
        ui->tabWidget->setCurrentIndex(0);
        QString msg;
        double dt;
        int IPA, IPV, IPD;
        IPA = peakLoc(eqs->getAcc(),eqs->getN());
        IPV = peakLoc(eqs->getVel(),eqs->getN());
        IPD = peakLoc(eqs->getDisp(),eqs->getN());
        dt = eqs->getDt();

        msg = QString("PA=%1@%4  PV=%2@%5  PD=%3@%6")
                .arg(eqs->getAcc()[IPA])
                .arg(eqs->getVel()[IPV])
                .arg(eqs->getDisp()[IPD])
                .arg(dt*(IPA)).arg(dt*(IPV)).arg(dt*(IPD));
        ui->statusBar->showMessage(msg);
    }


}

void MainWindow::plotRES()
{
    QCustomPlot *qplot = ui->ViewRES;
	qplot->plotLayout()->clear();
    qplot->clearGraphs();

    QCPAxisRect *plotAcc = new QCPAxisRect(qplot);
    QCPAxisRect *plotVel = new QCPAxisRect(qplot);
    QCPAxisRect *plotDsp = new QCPAxisRect(qplot);

    qplot->plotLayout()->addElement(0, 0, plotAcc);
    qplot->plotLayout()->addElement(1, 0, plotVel);
    qplot->plotLayout()->addElement(2, 0, plotDsp);

    QCPMarginGroup *marginGroup = new QCPMarginGroup(qplot);

    plotAcc->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);
    plotVel->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);
    plotDsp->setMarginGroup(QCP::msLeft | QCP::msRight, marginGroup);

    foreach(QCPAxisRect *rect, qplot->axisRects())
    {
        rect->axis(QCPAxis::atBottom)->setRange(0.0, 40.0);
        rect->axis(QCPAxis::atLeft)->setRange(-1.0, 1.0);
        foreach(QCPAxis *axis, rect->axes())
        {
            axis->setLayer("axes");
            axis->grid()->setLayer("grid");
        }
    }

    plotAcc->axis(QCPAxis::atLeft)->setLabel(tr("Acceleration"));
    plotVel->axis(QCPAxis::atLeft)->setLabel(tr("Velocity"));
    plotDsp->axis(QCPAxis::atLeft)->setLabel(tr("Displacement"));
    plotDsp->axis(QCPAxis::atBottom)->setLabel(tr("Time"));

    qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atLeft));
    qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atLeft));
    qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atLeft));

    qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atLeft));
    qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atLeft));
    qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atLeft));

    qplot->addGraph(plotAcc->axis(QCPAxis::atBottom), plotAcc->axis(QCPAxis::atLeft));
    qplot->addGraph(plotVel->axis(QCPAxis::atBottom), plotVel->axis(QCPAxis::atLeft));
    qplot->addGraph(plotDsp->axis(QCPAxis::atBottom), plotDsp->axis(QCPAxis::atLeft));

    QPen penRaw = QPen(Qt::gray);
    qplot->graph(0)->setPen(penRaw);
    qplot->graph(1)->setPen(penRaw);
    qplot->graph(2)->setPen(penRaw);

    QPen penTar = QPen(QColor(46, 139, 87));
    qplot->graph(6)->setPen(penTar);
    qplot->graph(7)->setPen(penTar);
    qplot->graph(8)->setPen(penTar);

    QVector<double> Time = eqs->qGetT();

    if (ui->INPUTON->isChecked())
    {
        QVector<double> Acc = eqs->qGetAcc();
        QVector<double> Vel = eqs->qGetVel();
        QVector<double> Dsp = eqs->qGetDisp();

        qplot->graph(0)->setData(Time, Acc);
        qplot->graph(1)->setData(Time, Vel);
        qplot->graph(2)->setData(Time, Dsp);

    }

    QVector<double> Ra = eqs->qGetRa();
    QVector<double> Rv = eqs->qGetRv();
    QVector<double> Rd = eqs->qGetRd();

    qplot->graph(3)->setData(Time, Ra);
    qplot->graph(4)->setData(Time, Rv);
    qplot->graph(5)->setData(Time, Rd);

	qplot->rescaleAxes();
	qplot->replot();
    ui->tabWidget->setCurrentWidget(ui->tabRes);

    int IPA = peakLoc(Ra.data(), Ra.count());
    int IPV = peakLoc(Rv.data(), Rv.count());
    int IPD = peakLoc(Rd.data(), Rd.count());
    double dt = eqs->getDt();

    QString msg;
    msg = QString("PA=%1@%4  PV=%2@%5  PD=%3@%6").arg(Ra[IPA]).arg(Rv[IPV]).arg(Rd[IPD])
        .arg(dt*(IPA)).arg(dt*(IPV)).arg(dt*(IPD));

    ui->statusBar->showMessage(msg);
}

void MainWindow::plotEnergy()
{
    double **e = eqs->getEnergy();
    int n = eqs->getN();

    QVector<double> t   = eqs->qGetT();
    QVector<double> Ek  = A2QV(e[0],n);
    QVector<double> Ez  = A2QV(e[1],n);
    QVector<double> Es  = A2QV(e[2],n);
    QVector<double> Eh  = A2QV(e[3],n);
    QVector<double> Ein = A2QV(e[4],n);

    QVector<double> Ez_Eh(n);
    QVector<double> Ez_Eh_Es(n);

    for (int i=0; i<n; ++i) {
        Ez_Eh[i] = Ez[i] + Eh[i];
        Ez_Eh_Es[i] = Ez_Eh[i] + Es[i];
    }

    QList< QVector<double> > E;
    E << Ez << Ez_Eh << Ez_Eh_Es << Ein;

    QCustomPlot *qplot = ui->ViewEnergy;
    qplot->clearGraphs();
    qplot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);

    qplot->legend->setVisible(true);

    QList<QCPGraph *> GR;
    GR << qplot->addGraph()
       << qplot->addGraph()
       << qplot->addGraph()
       << qplot->addGraph();

    QList<QColor> COLOR;
	COLOR << Qt::yellow << Qt::blue << Qt::red << Qt::darkGreen;

    QStringList NAME;
    NAME << tr("Viscous Damping Energy")
         << tr("Hysteretic Energy")
         << tr("Elastic Energy")
         << tr("Kinetic Energy");

    for (int i=0; i<4; ++i) {

        GR[i]->setName(NAME[i]);
        GR[i]->setPen(QPen(COLOR[i]));

        COLOR[i].setAlpha(100);
        GR[i]->setBrush(QBrush(COLOR[i]));

        GR[i]->setData(t,E[i]);
        if (i>0) {
            GR[i]->setChannelFillGraph(GR[i-1]);
            GR[i]->rescaleAxes(true);
        }
        else
            GR[i]->rescaleAxes();
    }

    qplot->replot();

//    ui->tabWidget->setCurrentWidget(ui->tabEnergy);
}

void MainWindow::plotHyst()
{
    QCustomPlot *qplot = ui->ViewHyst;
    qplot->clearPlottables();

//    qplot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
//    qplot->legend->setVisible(true);

    QVector<double> rd = eqs->qGetRd();
    QVector<double> rF = eqs->qGetRF();

    double xrl, xru, dx;
	double yrl, yru, dy;

    autoScale(rd.data(), eqs->getN(), xrl, xru, dx);
	autoScale(rF.data(), eqs->getN(), yrl, yru, dy);

    qplot->xAxis->setAutoTickStep(false);
    qplot->yAxis->setAutoTickStep(false);

    qplot->xAxis->setRange(xrl, xru);
    qplot->yAxis->setRange(yrl, yru);

    qplot->xAxis->setTickStep(dx);
    qplot->yAxis->setTickStep(dy);

    QCPCurve *cv = new QCPCurve(qplot->xAxis, qplot->yAxis);
    qplot->addPlottable(cv);

    cv->setData(rd,rF);

    qplot->replot();
//    ui->tabWidget->setCurrentWidget(ui->tabHyst);

}

void MainWindow::confirm()
{
	eqs->confirm();
	this->showDriftMsg();
    this->plotTH();
}

void MainWindow::saveTH()
{
    QString filename, sf;
    QStringList filters;
    QFileDialog fdialog(this);

    filters << "Text file (*.txt)" << "CSV file (*.csv)";

    fdialog.setAcceptMode(QFileDialog::AcceptSave);
    fdialog.setNameFilters(filters);
    fdialog.setWindowTitle(tr("Save"));
    fdialog.setDirectory(workDir);
    fdialog.selectFile(eqsName+"-TH");

    if (fdialog.exec()) {
        filename = fdialog.selectedFiles()[0];
        sf = fdialog.selectedNameFilter();
        if (sf == filters[0])
            eqs->savetxt(filename);
        else if (sf == filters[1])
            eqs->savecsv(filename);
    }
}

void MainWindow::setupSP()
{
    QString DR = ui->DR->text().trimmed();
    QStringList DRL = DR.split(",",QString::SkipEmptyParts);
    ui->CDR->clear();
    ui->CDR->addItems(DRL);
    int NSP = DRL.count();
    double *Zeta = new double[NSP];
    for (int i = 0; i < NSP; i++) Zeta[i] = DRL[i].toDouble();

    int np = ui->NP->value();
    double P1 = ui->TS->value();
    double P2 = ui->TL->value();
    int dm = ui->DM->currentIndex();
    int sm = ui->SM->currentIndex()+1;
    bool pseudo = ui->PD->isChecked();

    if (dm < 3) eqs->setupSP(NSP,np,P1,P2,dm,sm,pseudo,Zeta);
    else {
        double *p = pdw->getPeriods();
        eqs->setupSP(NSP,np,p,sm,pseudo,Zeta);
    }
}

void MainWindow::plotSPA()
{
    QCustomPlot *qplot = ui->ViewSPA;
    qplot->clearGraphs();

    qplot->legend->setVisible(true);

    if (ui->XS_LOG_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        qplot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignHCenter);
    }
    else if (ui->XS_LIN_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLinear);
        qplot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignRight);
    }

    for (int i = 0; i < eqs->getNsp(); i++)
    {
        QCPGraph *gr = qplot->addGraph();
        QPen pen = autoPen(i);
        gr->setPen(pen);
        Spectra *spi = eqs->getSP(i);

        gr->setData(spi->qGetP(),spi->qGetSPA());

        gr->setName(tr("Damping Ratio: %1%").arg((int)(spi->getZeta()*100.0)));
    }

    qplot->rescaleAxes();
	qplot->xAxis->setRangeLower(0.01);
	qplot->yAxis->setRangeLower(0.0);
    qplot->replot();

    ui->tabWidget->setCurrentIndex(2);

}

void MainWindow::plotSPAi()
{
    QCustomPlot *qplot = ui->ViewSPA;

    if (ui->XS_LOG_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        qplot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignHCenter);
    }
    else if (ui->XS_LIN_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLinear);
        qplot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignRight);
    }

    int i = ui->CDR->currentIndex();
    QCPGraph *gr_pre = qplot->addGraph();
    QCPGraph *gr_target = qplot->addGraph();

    QPen pen = QPen(Qt::gray);
    pen.setWidthF(1.5);
    gr_pre->setPen(pen);

    Spectra *spi = eqs->getSP(i);

    gr_pre->setData(spi->qGetP(),spi->qGetSPA());

    gr_pre->setName(tr("Before Fitting") + " (" + tr("Damping Ratio: %1%").arg((int)(spi->getZeta()*100.0)) + ")");

    pen.setStyle(Qt::DashLine);
    pen.setColor(Qt::blue);
    pen.setWidth(2);
    gr_target->setPen(pen);
    gr_target->setData(spi->qGetP(),spi->qGetSPT());
    gr_target->setName(tr("Target SPA") + " (" + tr("Damping Ratio: %1%").arg((int)(spi->getZeta()*100.0)) + ")");

    qplot->rescaleAxes();
    qplot->xAxis->setRangeLower(0.01);
    qplot->yAxis->setRangeLower(0.0);
    qplot->replot();

    ui->tabWidget->setCurrentWidget(ui->tabSPA);

}

void MainWindow::plotSPT()
{
    QCustomPlot *qplot = ui->ViewSPA;
    qplot->clearGraphs();

    if (ui->XS_LOG_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    }
    else if (ui->XS_LIN_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLinear);
    }

    bool rs = false;

    QCPGraph *gr;
    QPen pen;
    for (int i = 0; i < eqs->getNsp(); i++)
    {
        gr = qplot->addGraph();
        pen = autoPen(i);
        pen.setStyle(Qt::DashLine);
        gr->setPen(pen);
        Spectra *spi = eqs->getSP(i);

        gr->setData(spi->qGetP(),spi->qGetSPT());

        gr->rescaleAxes(rs);
        rs = rs + true;

        gr->setName(tr("Damping Ratio: %1%").arg((int)(spi->getZeta()*100.0)) + " (" + tr("Target SPA") + ")");

        if (ui->SPAON->isChecked()) {
            QCPGraph *gra = qplot->addGraph();
            pen.setStyle(Qt::SolidLine);
            gra->setPen(pen);
            gra->setData(spi->qGetP(),spi->qGetSPA());
            gra->rescaleAxes(rs);

            gra->setName(tr("Damping Ratio: %1%").arg((int)(spi->getZeta()*100.0)));
        }

    }
    qplot->xAxis->setRangeLower(0.01);
    qplot->yAxis->setRangeLower(0.0);
    qplot->replot();

	ui->tabWidget->setCurrentWidget(ui->tabSPA);

}

void MainWindow::fillSPTable(int k)
{
    Spectra *spi = eqs->getSP(k);

    double *P = spi->getP();
    double *SPA = spi->getSPA();
    double *SPV = spi->getSPV();
    double *SPD = spi->getSPD();
    double *SPE = spi->getSPE();
    double *SPT = spi->getSPT();

    dataTableSP->clearContents();
    dataTableSP->setRowCount(spi->getNP());

    for(int i=0; i<spi->getNP(); ++i) {
        dataTableSP->setItem(i,0,new QTableWidgetItem(QString::number(P[i])));
        dataTableSP->setItem(i,1,new QTableWidgetItem(QString::number(SPA[i])));
        dataTableSP->setItem(i,2,new QTableWidgetItem(QString::number(SPV[i])));
        dataTableSP->setItem(i,3,new QTableWidgetItem(QString::number(SPD[i])));
        dataTableSP->setItem(i,4,new QTableWidgetItem(QString::number(SPE[i])));
        dataTableSP->setItem(i,5,new QTableWidgetItem(QString::number(SPT[i])));
    }
}

void MainWindow::plotSP()
{
    QCustomPlot *qplot = ui->ViewSP;

    qplot->clearGraphs();
    qplot->plotLayout()->clear();

    QCPAxisRect *axisSPA = new QCPAxisRect(qplot);
    QCPAxisRect *axisSPV = new QCPAxisRect(qplot);
    QCPAxisRect *axisSPD = new QCPAxisRect(qplot);
    QCPAxisRect *axisSPE = new QCPAxisRect(qplot);

    qplot->plotLayout()->addElement(0, 0, axisSPA);
    qplot->plotLayout()->addElement(0, 1, axisSPV);
    qplot->plotLayout()->addElement(1, 0, axisSPD);
    qplot->plotLayout()->addElement(1, 1, axisSPE);

    QCPMarginGroup *marginGroup = new QCPMarginGroup(qplot);

    axisSPA->setMarginGroup(QCP::msAll, marginGroup);
    axisSPV->setMarginGroup(QCP::msAll, marginGroup);
    axisSPD->setMarginGroup(QCP::msAll, marginGroup);
    axisSPE->setMarginGroup(QCP::msAll, marginGroup);

    QCPLegend *legendSPA = new QCPLegend;
    QCPLegend *legendSPV = new QCPLegend;
    QCPLegend *legendSPD = new QCPLegend;
    QCPLegend *legendSPE = new QCPLegend;

    QBrush brush = QBrush(QColor(255,255,255,150));

    legendSPA->setBrush(brush);
    legendSPV->setBrush(brush);
    legendSPD->setBrush(brush);
    legendSPE->setBrush(brush);

    if (ui->XS_LOG_SP->isChecked()) {
        axisSPA->insetLayout()->addElement(legendSPA, Qt::AlignBottom|Qt::AlignHCenter);
        axisSPV->insetLayout()->addElement(legendSPV, Qt::AlignTop|Qt::AlignLeft);
        axisSPD->insetLayout()->addElement(legendSPD, Qt::AlignTop|Qt::AlignLeft);
        axisSPE->insetLayout()->addElement(legendSPE, Qt::AlignTop|Qt::AlignLeft);
    }
    else {
        axisSPA->insetLayout()->addElement(legendSPA, Qt::AlignTop|Qt::AlignRight);
        axisSPV->insetLayout()->addElement(legendSPV, Qt::AlignBottom|Qt::AlignRight);
        axisSPD->insetLayout()->addElement(legendSPD, Qt::AlignTop|Qt::AlignLeft);
        axisSPE->insetLayout()->addElement(legendSPE, Qt::AlignBottom|Qt::AlignRight);
    }

    legendSPA->setLayer("legend");
    legendSPV->setLayer("legend");
    legendSPD->setLayer("legend");
    legendSPE->setLayer("legend");

    legendSPA->setBorderPen(QPen(Qt::transparent));
    legendSPV->setBorderPen(QPen(Qt::transparent));
    legendSPD->setBorderPen(QPen(Qt::transparent));
    legendSPE->setBorderPen(QPen(Qt::transparent));

    qplot->setAutoAddPlottableToLegend(false);

    foreach(QCPAxisRect *rect, qplot->axisRects())
    {
        rect->axis(QCPAxis::atBottom)->setRange(0.01, 10.0);
        if (ui->XS_LOG_SP->isChecked())
        {
            rect->axis(QCPAxis::atBottom)->setScaleType(QCPAxis::stLogarithmic);
        }
        else if (ui->XS_LIN_SP->isChecked())
        {
            rect->axis(QCPAxis::atBottom)->setScaleType(QCPAxis::stLinear);
        }
        rect->axis(QCPAxis::atBottom)->grid()->setSubGridVisible(true);
        rect->axis(QCPAxis::atBottom)->setLabel(tr("Period"));
        rect->axis(QCPAxis::atLeft)->setRange(0.0, 2.5);

        foreach(QCPAxis *axis, rect->axes())
        {
            axis->setLayer("axes");
            axis->grid()->setLayer("grid");
        }
    }

    axisSPA->axis(QCPAxis::atLeft)->setLabel(tr("SPA"));
    axisSPV->axis(QCPAxis::atLeft)->setLabel(tr("SPV"));
    axisSPD->axis(QCPAxis::atLeft)->setLabel(tr("SPD"));
    axisSPE->axis(QCPAxis::atLeft)->setLabel(tr("SPE"));

    bool rs = false;

    for (int i = 0; i < eqs->getNsp(); ++i)
    {
        QCPGraph *grSPA = qplot->addGraph(axisSPA->axis(QCPAxis::atBottom), axisSPA->axis(QCPAxis::atLeft));
        QCPGraph *grSPV = qplot->addGraph(axisSPV->axis(QCPAxis::atBottom), axisSPV->axis(QCPAxis::atLeft));
        QCPGraph *grSPD = qplot->addGraph(axisSPD->axis(QCPAxis::atBottom), axisSPD->axis(QCPAxis::atLeft));
        QCPGraph *grSPE = qplot->addGraph(axisSPE->axis(QCPAxis::atBottom), axisSPE->axis(QCPAxis::atLeft));
        
        grSPA->setPen(autoPen(i));
        grSPV->setPen(autoPen(i));
        grSPD->setPen(autoPen(i));
        grSPE->setPen(autoPen(i));

        Spectra *spi = eqs->getSP(i);

        QVector<double> P = spi->qGetP();
        QVector<double> SPA = spi->qGetSPA();
        QVector<double> SPV = spi->qGetSPV();
        QVector<double> SPD = spi->qGetSPD();
        QVector<double> SPE = spi->qGetSPE();

        grSPA->setData(P,SPA);
        grSPV->setData(P,SPV);
        grSPD->setData(P,SPD);
        grSPE->setData(P,SPE);

        if (i>0) rs = true;
        grSPA->rescaleValueAxis(rs);
        grSPV->rescaleValueAxis(rs);
        grSPD->rescaleValueAxis(rs);
        grSPE->rescaleValueAxis(rs);

        QString name = tr("Damping Ratio: %1%").arg((int)(spi->getZeta()*100.0));
        grSPA->setName(name);
        grSPV->setName(name);
        grSPD->setName(name);
        grSPE->setName(name);

        legendSPA->addItem(new QCPPlottableLegendItem(legendSPA, grSPA));
        legendSPV->addItem(new QCPPlottableLegendItem(legendSPV, grSPV));
        legendSPD->addItem(new QCPPlottableLegendItem(legendSPD, grSPD));
        legendSPE->addItem(new QCPPlottableLegendItem(legendSPE, grSPE));
    }

    axisSPA->axis(QCPAxis::atLeft)->setRangeLower(0.0);
    axisSPV->axis(QCPAxis::atLeft)->setRangeLower(0.0);
    axisSPD->axis(QCPAxis::atLeft)->setRangeLower(0.0);
    axisSPE->axis(QCPAxis::atLeft)->setRangeLower(0.0);

    qplot->replot();

    ui->tabWidget->setCurrentIndex(3);

}

void MainWindow::plotFFT()
{
    QCustomPlot *qplot = ui->ViewFFT;
    qplot->clearGraphs();
    qplot->clearPlottables();

    ui->Amp->setEnabled(true);
    ui->Ang->setEnabled(true);
    ui->DAng->setEnabled(true);

    QCPGraph *gr;
    if (ui->Amp->isChecked()) {
        gr = qplot->addGraph();
        gr->setData(eqs->qGetFreqs(),eqs->qGetAmpf());
        gr->setLineStyle(QCPGraph::lsImpulse);
        qplot->xAxis->setRangeReversed(false);
        qplot->xAxis->setAutoTicks(true);
        qplot->xAxis->setAutoTickLabels(true);
        qplot->xAxis->setAutoTickStep(true);
        qplot->yAxis->setAutoTickStep(true);
        qplot->xAxis->setAutoSubTicks(true);
        qplot->yAxis->setAutoSubTicks(true);
        qplot->xAxis->setLabel(tr("Frequency (Hz)"));
        qplot->yAxis->setLabel(tr("Fourier Amplitude"));
    }
    else if (ui->Ang->isChecked()) {

        QVector<double> dAng = eqs->qGetAngf();
        QCPBars *bars = new QCPBars(qplot->xAxis,qplot->yAxis);

        int N = 100;
        QVector<double> x(N), y(N);

        hist(dAng,0.0,2.0*PI,N,x,y);
        bars->setData(x,y);
        bars->setWidth(2.0*PI/N);
        qplot->addPlottable(bars);

        qplot->xAxis->setRangeReversed(false);
        qplot->xAxis->setAutoTicks(false);
        qplot->xAxis->setAutoTickLabels(false);

        QVector<double> phi;
        phi << 0.0 << PI << PI*2;
        QVector<QString> label;
        label << "0" << QString::fromUtf8("π")
                     << QString::fromUtf8("2π");
        qplot->xAxis->setTickVector(phi);
        qplot->xAxis->setTickVectorLabels(label);
        qplot->xAxis->setAutoSubTicks(false);
        qplot->xAxis->setSubTickCount(0);
        qplot->xAxis->setLabel(tr("Phase Angle"));

        qplot->yAxis->setAutoTickStep(false);
        int tstp = 10*(int)(max(y.data(),N)/50);
        tstp = tstp<10? 10:tstp;
        qplot->yAxis->setTickStep(tstp);
        qplot->yAxis->setAutoSubTicks(false);
        qplot->yAxis->setSubTickCount(1);
        qplot->yAxis->setLabel(tr("Count"));
    }
    else if (ui->DAng->isChecked()) {
        QVector<double> dAng = eqs->qGetDAngf();
        QCPBars *bars = new QCPBars(qplot->xAxis,qplot->yAxis);

        int N = 100;
        QVector<double> x(N), y(N);

        hist(dAng,-2.0*PI,0.0,N,x,y);
        bars->setData(x,y);
        bars->setWidth(2.0*PI/N);
        qplot->addPlottable(bars);
        qplot->xAxis->setRangeReversed(true);
        qplot->xAxis->setAutoTicks(false);
        qplot->xAxis->setAutoTickLabels(false);

        QVector<double> phi;
        phi << 0.0 << -PI << -PI*2;
        QVector<QString> label;
        label << "0" << QString::fromUtf8("-π")
                     << QString::fromUtf8("-2π");
        qplot->xAxis->setTickVector(phi);
        qplot->xAxis->setTickVectorLabels(label);
        qplot->xAxis->setAutoSubTicks(false);
        qplot->xAxis->setSubTickCount(1);
        qplot->xAxis->setLabel(tr("Phase Difference"));

        qplot->yAxis->setAutoTickStep(false);
        int tstp = 10*(int)(max(y.data(),N)/50);
        tstp = tstp<10? 10:tstp;
        qplot->yAxis->setTickStep(tstp);
        qplot->yAxis->setAutoSubTicks(false);
        qplot->yAxis->setSubTickCount(1);
        qplot->yAxis->setLabel(tr("Count"));

    }



    qplot->xAxis->setScaleType(QCPAxis::stLinear);

    ui->XS_LOG_SP->setChecked(false);
    ui->XS_LIN_SP->setChecked(true);

//    if (ui->XS_LOG_SP->isChecked())
//    {
//        qplot->xAxis->setScaleType(QCPAxis::stLogarithmic);
//    }
//    else if (ui->XS_LIN_SP->isChecked())
//    {
//        qplot->xAxis->setScaleType(QCPAxis::stLinear);
//    }

    qplot->rescaleAxes();

    qplot->replot();

    ui->tabWidget->setCurrentIndex(1);

}

void MainWindow::plotPSD()
{
    QCustomPlot *qplot = ui->ViewFFT;
    qplot->clearGraphs();

    qplot->yAxis->setLabel(tr("PSD"));

    QCPGraph *gr = qplot->addGraph();
    gr->setData(eqs->qGetFpsd(),eqs->qGetPsd());

//    gr->setLineStyle(QCPGraph::lsImpulse);
//    gr->setPen(QPen(Qt::blue,1.0));

    if (ui->XS_LOG_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    }
    else if (ui->XS_LIN_SP->isChecked())
    {
        qplot->xAxis->setScaleType(QCPAxis::stLinear);
    }

    gr->rescaleAxes();

    qplot->replot();

    ui->Amp->setDisabled(true);
    ui->Ang->setDisabled(true);

    ui->tabWidget->setCurrentIndex(1);
}

void MainWindow::SPXScaleChanged()
{

    switch (ui->tabWidget->currentIndex()) {
    case 1:
        if (ui->XS_LOG_SP->isChecked())
        {
            ui->ViewFFT->xAxis->setScaleType(QCPAxis::stLogarithmic);
        }
        else if (ui->XS_LIN_SP->isChecked())
        {
            ui->ViewFFT->xAxis->setScaleType(QCPAxis::stLinear);
        }
        ui->ViewFFT->replot();
        break;
    case 2:
        if (ui->XS_LOG_SP->isChecked())
        {
            ui->ViewSPA->xAxis->setScaleType(QCPAxis::stLogarithmic);
            ui->ViewSPA->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignHCenter);
        }
        else if (ui->XS_LIN_SP->isChecked())
        {
            ui->ViewSPA->xAxis->setScaleType(QCPAxis::stLinear);
            ui->ViewSPA->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignRight);
        }
        ui->ViewSPA->replot();
        break;
    case 3:
        this->plotSP();
    default:
        break;
    }
}

/*
void MainWindow::localAdjust()
{
    int ind1 = ui->Ind1->value();
    int ind2 = ui->Ind2->value();
    int nl = ind2 - ind1 + 1;
    double *a = new double[nl];
    double *eqsa = eqs->getAcc();
    for (int i = 0; i < nl; i++) a[i] = eqsa[i + ind1];
    EQSignal *leqs = new EQSignal(a, nl, eqs->getDt(), *(eqs->getVel() + ind1), *(eqs->getDisp() + ind1));
    bool raw = true;
    if (ui->checkFilt->isChecked())
    {
        double f1 = ui->Freq1->value();
        double f2 = ui->Freq2->value();
        int ftype = ui->FilterType->currentIndex();
        int order = ui->FOrder->value();
        leqs->filt(ftype, order, f1, f2);
        raw = false;
    }
    if (ui->checkDetrend->isChecked())
    {
        int oh = ui->Order2->value();
        int ol = ui->Order1->value();
        ol = max(ol, 1);
        int method = ui->DetrendMethod->currentIndex();
        leqs->detrend(method, oh, ol, raw);
        raw = false;
    }
    if (ui->checkAlign->isChecked())
    {
        int oh = ui->TDispTrend->value();
        int ol = 0;
        int ntp = ui->NTP->value();
        int method = ui->AlignMethod->currentIndex();
        bool EWZ = ui->EndWZ->isChecked();
        leqs->align(method, ntp, oh, ol, raw, EWZ);
    }

    double *leqsa = leqs->getAcc();

    for (int i = 0; i < nl; i++) eqsa[i + ind1] = leqsa[i];

    this->showDriftMsg();
    eqs->a2vd();
    this->plotTH();

    delete[] a;
}
*/

void MainWindow::showDriftMsg()
{
	double dr0 = eqs0->getDR();
	double ar0 = eqs0->getAR();
	double dr = eqs->getDR();
	double ar = eqs->getAR();

	QString msg;
	msg.sprintf("DR: %.4g>>%.4g\t\tAR: %.4g>>%.4g", dr0, dr, ar0, ar);
	ui->statusBar->showMessage(msg);
}

void MainWindow::on_Norm_clicked()
{
    eqs->norm();
    eqs0->norm();
    this->plotTH();
}

void MainWindow::on_Reload_clicked()
{
    QString filename = ui->URL->text();

    if (filename.endsWith(".txt"))
    {
        double DT;
        DT = QInputDialog::getDouble(this, tr("Time Interval"), "dt = ", 0.02, 0.0, 10.0, 3);
        this->readtxt(filename, DT);
    }
    else if (filename.endsWith(".at2"))
    {
        this->readnga(filename);
    }

}

void MainWindow::on_Filt_clicked()
{
    double f1 = ui->Freq1->value();
    double f2 = ui->Freq2->value();
    int ftype = ui->FilterType->currentIndex();
    int order = ui->FOrder->value();

    eqs->filt(ftype, order, f1, f2);
    this->showDriftMsg();
    this->plotTH();
}

void MainWindow::on_Detrend_clicked()
{
    int oh = ui->Order2->value();
    int ol = ui->Order1->value();
    int method = ui->DetrendMethod->currentIndex();

    eqs->detrend(method, oh, ol);
    this->showDriftMsg();

    this->plotTH();
}

void MainWindow::on_Align_clicked()
{
    int oh = ui->TDispTrend->value();
    int ol = 0;
    int ntp = ui->NTP->value();
    int method = ui->AlignMethod->currentIndex();
    bool EWZ = ui->EndWZ->isChecked();
    eqs->align(method, ntp, oh, ol, true, EWZ);
    this->showDriftMsg();

    this->plotTH();

    if (ui->TargetON->isChecked())
    {
        QCustomPlot *qplot = ui->ViewTH;
        QVector<double> T = eqs->qGetT();
        QVector<double> TA = eqs0->qGetTa();
        QVector<double> TV = eqs0->qGetTv();
        QVector<double> TD = eqs0->qGetTd();

        qplot->graph(6)->setData(T, eqs->qGetTa());
        qplot->graph(6)->rescaleValueAxis(true);
        qplot->graph(7)->setData(T, eqs->qGetTv());
        qplot->graph(7)->rescaleValueAxis(true);
        qplot->graph(8)->setData(T, eqs->qGetTd());
        qplot->graph(8)->rescaleValueAxis(true);
        
        qplot->replot();
    }
}

void MainWindow::on_Adjust_clicked()
{
    bool raw = true;
    if (ui->checkFilt->isChecked())
    {
        double f1 = ui->Freq1->value();
        double f2 = ui->Freq2->value();
        int ftype = ui->FilterType->currentIndex();
        int order = ui->FOrder->value();
        eqs->filt(ftype, order, f1, f2);
        raw = false;
    }
    if (ui->checkDetrend->isChecked())
    {
        int oh = ui->Order2->value();
        int ol = ui->Order1->value();
        int method = ui->DetrendMethod->currentIndex();
        eqs->detrend(method, oh, ol, raw);
        raw = false;
    }
    if (ui->checkAlign->isChecked())
    {
        int oh = ui->TDispTrend->value();
        int ol = 0;
        int ntp = ui->NTP->value();
        int method = ui->AlignMethod->currentIndex();
        bool EWZ = ui->EndWZ->isChecked();
        eqs->align(method, ntp, oh, ol, raw, EWZ);
    }
    this->showDriftMsg();
    this->plotTH();
}

void MainWindow::on_SetTrimEdge_clicked()
{
    int method = ui->STEMethod->currentIndex();
    double thd1 = ui->Threshold1->value();
    double thd2 = ui->Threshold2->value();
    bool EZ = ui->EZ->isChecked();

    int *edges = eqs->autoTrimEdges(method,thd1,thd2,EZ);
    ui->Ind1->setValue(edges[0]);
    ui->Ind2->setValue(edges[1]);
}

void MainWindow::on_Trim_clicked()
{
    int ind1 = ui->Ind1->value();
    int ind2 = ui->Ind2->value();
    eqs->trim(ind1,ind2);
    eqs0->trim(ind1,ind2);

    this->plotTH();
}

void MainWindow::on_CalcFFT_clicked()
{
    eqs->calcFFT();
    this->plotFFT();
}

void MainWindow::on_CalcPSD_clicked()
{
    eqs->calcPSD(ui->OLR->value(),ui->WW->isChecked());
    this->plotPSD();
}

void MainWindow::on_CalcSPA_clicked()
{
    setupSP();
    eqs->calcSP(false);
    plotSPA();
    FIT_SPA = false;
}

void MainWindow::on_CalcSP_clicked()
{
    setupSP();
    eqs->calcSP(true);
    plotSP();
    FIT_SPA = false;
}

//void MainWindow::on_DR_textChanged(const QString &arg1)
//{

//    QString DR;
//    DR.append(arg1.trimmed());
//    QStringList DRL = DR.split(",");
//    ui->CDR->clear();
//    ui->CDR->addItems(DRL);

//    dr->clear();
//    dr->addItems(DRL);
//}

//void MainWindow::on_tabWidget_tabBarClicked(int index)
//{
//    switch (index) {
//    case 0:
//        ui->Paras->setCurrentIndex(0);
//        break;
//    case 1:
//        ui->Paras->setCurrentIndex(3);
//        break;
//    case 2:
//        ui->Paras->setCurrentIndex(4);
//        break;
//    case 3:
//        ui->Paras->setCurrentIndex(4);
//        break;
//    default:
//        break;
//    }
//}

void MainWindow::on_actionOpen_triggered()
{
    QString fileName;
    QStringList filters;
    QFileDialog fdialog(this);

    filters << "Text file (*.txt *.TXT)" << "NGA record file (*.at2 *.AT2)";

    fdialog.setNameFilters(filters);
    fdialog.setWindowTitle(tr("Open"));
    fdialog.setDirectory(workDir);

    if (fdialog.exec()) fileName = fdialog.selectedFiles()[0];

    if (fileName.isEmpty()) return;

    bool isTxt = fileName.endsWith(".txt") || (fileName.endsWith(".TXT"));
    bool isAt2 = fileName.endsWith(".at2") || (fileName.endsWith(".AT2"));

    QFileInfo fi(fileName);
    eqsName = fi.baseName();
    this->setWindowTitle("EQSignal -- " + eqsName);

    if (isTxt)
    {
        ui->URL->setText(fileName);
        double DT = QInputDialog::getDouble(this, eqsName + "'s " + tr("Time Interval"), "dt = ", 0.02, 0.0, 10.0, 3);
        this->readtxt(fileName, DT);

    }
    else if (isAt2)
    {
        ui->URL->setText(fileName);
        this->readnga(fileName);
    }

    lastFile = fileName;

}

void MainWindow::on_actionSaveData_triggered()
{
    switch (ui->tabWidget->currentIndex()) {
    case 0:
        this->saveTH();
        break;
//    case 1:
//        this->saveFFT();
//        break;
    case 2:
        this->saveSP();
        break;
    case 3:
        this->saveSP();
        break;
    default:
        break;
    }
}

void MainWindow::on_actionSaveFig_triggered()
{
    QString figName,autoName;
    autoName.append(workDir.path());
	autoName.append("/");
    QCustomPlot *qplot;

    switch (ui->tabWidget->currentIndex()) {
    case 0:
        qplot = ui->ViewTH;
        autoName.append(eqsName);
        autoName.append("-TH");
        break;
    case 1:
        qplot = ui->ViewFFT;
        autoName.append(eqsName);
        autoName.append("-FFT");
        break;
    case 2:
        qplot = ui->ViewSPA;
        autoName.append(eqsName);
        autoName.append("-SPA");
        break;
    case 3:
        qplot = ui->ViewSP;
        autoName.append(eqsName);
        autoName.append("-SP");
        break;
	case 4:
		qplot = ui->ViewRES;
		autoName.append(eqsName);
		autoName.append("-RES");
		break;
	case 5:
		qplot = ui->ViewEnergy;
		autoName.append(eqsName);
		autoName.append("-RES-Energy");
		break;
	case 6:
		qplot = ui->ViewHyst;
		autoName.append(eqsName);
		autoName.append("-RES-Hyst");
		break;
    default:
        qplot = ui->ViewTH;
        autoName.append(eqsName);
        autoName.append("-TH");
        break;
    }

	autoName.append(".png");
    figName = QFileDialog::getSaveFileName(this, tr("Save"), autoName, "*.png");

    qplot->setBackground(Qt::transparent);
    qplot->savePng(figName,800,600);
	qplot->setBackground(Qt::white);

//    QString filename, sf;
//    QStringList filters;
//    QFileDialog fdialog(this);

//    filters << "Portable Network Graphic Format (*.png)"
//            << "JPEG Format (*.jpg)"
//            << "Portable Document Format (*.pdf)";

//    fdialog.setAcceptMode(QFileDialog::AcceptSave);
//    fdialog.setNameFilters(filters);
//    fdialog.setWindowTitle("Save Figure");
//    fdialog.setDirectory(workDir);
//    fdialog.selectFile(eqsName+"-Acc");
}

void MainWindow::on_actionSaveAcc_triggered()
{
    QString filename, sf;
    QStringList filters;
    QFileDialog fdialog(this);

    double *time = eqs->getT();
    double *acc = eqs->getAcc();

    filters << "Text file (*.txt)" << "CSV file (*.csv)";

    fdialog.setAcceptMode(QFileDialog::AcceptSave);
    fdialog.setNameFilters(filters);
    fdialog.setWindowTitle(tr("Save"));
    fdialog.setDirectory(workDir);
    fdialog.selectFile(eqsName+"-Acc");

    if (fdialog.exec()) {
        filename = fdialog.selectedFiles()[0];
        sf = fdialog.selectedNameFilter();
        if (sf == filters[0]) {
            ofstream out(filename.toUtf8().data(), ios::out);
            for (int i = 0; i<eqs->getN(); ++i)
                if (saveAccWithTime)
                    out << time[i] << "\t" << acc[i] << endl;
                else
                    out << acc[i] << endl;
            out.close();
        }
        else if (sf == filters[1]) {
            ofstream out(filename.toUtf8().data(), ios::out);
            for (int i = 0; i<eqs->getN(); ++i)
                if (saveAccWithTime)
                    out << time[i] << ", " << acc[i] << endl;
                else
                    out << acc[i] << endl;
            out.close();
        }
    }
}

void MainWindow::fillTHTable()
{
    double *t = eqs->getT();
    double *a = eqs->getAcc();
    double *v = eqs->getVel();
    double *d = eqs->getDisp();

    dataTableTH->clearContents();
    dataTableTH->setRowCount(eqs->getN());

    for(int i=0; i<eqs->getN(); ++i) {
        dataTableTH->setItem(i,0,new QTableWidgetItem(QString::number(t[i])));
        dataTableTH->setItem(i,1,new QTableWidgetItem(QString::number(a[i])));
        dataTableTH->setItem(i,2,new QTableWidgetItem(QString::number(v[i])));
        dataTableTH->setItem(i,3,new QTableWidgetItem(QString::number(d[i])));
    }

//    dataTableTH->show();
}

void MainWindow::fillRESTable()
{
    double *t = eqs->getT();
    double *ra = eqs->getRa();
    double *rv = eqs->getRv();
    double *rd = eqs->getRd();

    dataTableRES->clearContents();
    dataTableRES->setRowCount(eqs->getN());

    for(int i=0; i<eqs->getN(); ++i) {
        dataTableRES->setItem(i,0,new QTableWidgetItem(QString::number(t[i])));
        dataTableRES->setItem(i,1,new QTableWidgetItem(QString::number(ra[i])));
        dataTableRES->setItem(i,2,new QTableWidgetItem(QString::number(rv[i])));
        dataTableRES->setItem(i,3,new QTableWidgetItem(QString::number(rd[i])));
    }

}

void MainWindow::genWave()
{
	
	int N = gwd->N->value();
    double dt = gwd->dt->value();
    double A0 = gwd->A0->value();
    double A = gwd->A->value();
    double T = gwd->T->value();
    double phi = gwd->phi->value();

    double *t = linspace(0.0,dt*N-dt,N);
    double *a = zeros(N);
	
	switch (gwd->wavetype->currentIndex())
	{
	case(0) :
		for (int i = 0; i < N; i++)
		{
			a[i] = A0 + A*sin(PI2 / T*t[i] + phi);
		}
		break;
	default:
		break;
	}

	eqs = new EQSignal(a,N,dt);
	eqs->a2vd();

	eqs0 = new EQSignal(a, N, dt);
	eqs0->a2vd();

	plotTH();

}

void MainWindow::showSPAErrorMsg()
{
    int i = ui->CDR->currentIndex();
    double Emax, Emean, CV;
    eqs->getSP(i)->fitError(Emax, Emean, CV);

    QString ErrInfo = tr("Before Fitting")
            + tr("Max Error: %1%, Mean Error: %2%")
            .arg(Emax*100.0).arg(Emean*100.0)
            + ", " + tr("CV: %1").arg(CV);
    ui->statusBar->showMessage(ErrInfo);
}

void MainWindow::on_GenSPT_clicked()
{
    if (ui->TSPT->currentIndex() == 0) {
        double Tg = ui->Tg->value();
        double PAF = ui->PAF->value();
        double SF = ui->SF->value();

        if (SF<0.0) SF = fabs(eqs->getPeakAcc());

        eqs->setSPT(Tg,PAF,SF);
    }
    this->plotSPT();
    FIT_SPA = false;
}

void MainWindow::on_actionViewData_triggered()
{
    int tabi = ui->tabWidget->currentIndex();

    if (tabi == 0) {
        this->fillTHTable();
        wth->show();
    }
    else if (tabi == 2 || tabi == 3) {
        this->fillSPTable();
        wsp->show();
    }
    else if (tabi == 4) {
        this->fillRESTable();
        wres->show();
    }

}

void MainWindow::on_DM_activated(int index)
{
    if (index == 3) {
        pdw->show();
        ui->NP->setEnabled(false);
        ui->TS->setEnabled(false);
        ui->TL->setEnabled(false);
    }
    else {
        ui->NP->setEnabled(true);
        ui->TS->setEnabled(true);
        ui->TL->setEnabled(true);
    }
}

void MainWindow::on_TSPT_activated(int index)
{
    switch (index) {
    case 1:
    case 2:
        sptw->show();
        ui->Tg->setEnabled(false);
        ui->PAF->setEnabled(false);
        ui->SF->setEnabled(false);
        break;
    default:
        ui->Tg->setEnabled(true);
        ui->PAF->setEnabled(true);
        ui->SF->setEnabled(true);
        break;
    }
}

void MainWindow::pdw_applyButton_clicked()
{
    double *p = pdw->getPeriods();
    int np = pdw->dataTable->rowCount();

    for (int i=0; i<eqs->getNsp(); ++i) eqs->getSP(i)->setP(p,np);

    ui->TS->setValue(p[0]);
    ui->TL->setValue(p[np-1]);
    ui->NP->setValue(np);
}

void MainWindow::sptw_applyButton_clicked()
{
    double *p = sptw->getPeriods();
    double *spt = sptw->getSPT();
    int np = sptw->dataTable->rowCount();

    for (int i=0; i<eqs->getNsp(); ++i) eqs->getSP(i)->setSPT(p,spt,np);

    ui->TS->setValue(p[0]);
    ui->TL->setValue(p[np-1]);
    ui->NP->setValue(np);
}

void MainWindow::on_SPFit_clicked()
{
    int i = ui->CDR->currentIndex();
    double tol = ui->Tol->value();
	double peak0 = ui->SF->value();
    int mit = ui->MIT->value();
    int fm = ui->FitMethod->currentIndex();
    int iter = 0;

    FIT_SPA = true;

    ui->ViewSPA->clearGraphs();

    if (ui->tabWidget->currentIndex() != 2)
        ui->tabWidget->setCurrentIndex(2);

    this->plotSPAi();

    Spectra *spi = eqs->getSP(i);
    double Emax, Emean, CV;
    double Emaxp, Emeanp, CVp;
    QString ErrInfo;
	QString msg;
    spi->fitError(Emax, Emean, CV);
    ErrInfo = tr("Before Fitting")
            + tr("Max Error: %1%, Mean Error: %2%")
            .arg(Emax*100.0).arg(Emean*100.0)
            + ", " + tr("CV: %1").arg(CV);
    ui->statusBar->showMessage(ErrInfo);

	Emaxp = Emax;
	Emeanp = Emean;
    CVp = CV;

    if (Emax <= 1.2*tol) {
        msg = tr("Error is less than Tolerance. Need no Fitting!");
		QMessageBox::information(0, tr("EQSignal"), msg);
		return;
	}
    msg = tr("Fitting, please waiting ...");
    ui->statusBar->showMessage(msg);

    QCustomPlot *qplot = ui->ViewSPA;
    QCPGraph *gr_post = qplot->addGraph();
    QPen pen(Qt::blue);
    pen.setWidthF(1.5);
    gr_post->setPen(pen);
    gr_post->setName(tr("After Fitting") + " (" + tr("Damping Ratio: %1%").arg((int)(spi->getZeta()*100.0)) + ")");

    clock_t ct = clock();
	if (fm == 0)
	{
		eqs->fitSP(i, tol, mit, fm, peak0);
		eqs->calcSP(i);
        spi->fitError(Emax, Emean, CV);

        gr_post->setData(spi->qGetP(), spi->qGetSPA());
        qplot->rescaleAxes();
        qplot->xAxis->setRangeLower(0.01);
        qplot->yAxis->setRangeLower(0.0);
        qplot->replot();

        if (Emax <= 1.2*tol) {
            msg = tr("Spectrum Fitting Converged not more than %1 iterations!").arg(mit)
                + tr("\nTotal Consumed Time: %1 s").arg((double)(clock()-ct)/CLOCKS_PER_SEC);
            QMessageBox::information(0, tr("EQSignal"), msg);
        }
        else {
            msg = tr("Spectrum Fitting not Converged After %1 iterations!").arg(mit)
                    + tr("\nTotal Consumed Time: %1 s").arg((double)(clock()-ct)/CLOCKS_PER_SEC);
            QMessageBox::warning(0, tr("EQSignal"), msg);
        }

	}
    else if (fm == 1)
    {
        ui->progressBar->setMaximum(mit);
        ui->progressBar->show();
		
        while (Emax > 1.2*tol && iter<mit)
		{
            iter ++;
            eqs->fitSP(i, tol, 1, fm, peak0);
            eqs->calcSP(i);
            spi->fitError(Emax, Emean, CV);

            gr_post->setData(spi->qGetP(), spi->qGetSPA());            
            qplot->rescaleAxes();
            qplot->xAxis->setRangeLower(0.01);
            qplot->yAxis->setRangeLower(0.0);
            qplot->replot();

            ui->progressBar->setValue(iter);
		}

        if (ui->progressBar->value()<mit)
            ui->progressBar->setValue(mit);
        ui->progressBar->hide();

        if (Emax <= 1.2*tol) {
            msg = tr("Spectrum Fitting Converged After %1 iterations!").arg(iter)
                    + tr("\nTotal Consumed Time: %1 s").arg((double)(clock()-ct)/CLOCKS_PER_SEC);
            QMessageBox::information(0, tr("EQSignal"), msg);
        }
        else {
            msg = tr("Spectrum Fitting not Converged After %1 iterations!").arg(iter)
                    + tr("\nTotal Consumed Time: %1 s").arg((double)(clock()-ct)/CLOCKS_PER_SEC);
            QMessageBox::warning(0, tr("EQSignal"), msg);
        }

	}
    else
    {
        ui->progressBar->setMaximum(mit);
        ui->progressBar->show();
        eqs->fitSP(i, tol, 5, 0, peak0);

        while (Emax > 1.2*tol && iter<mit)
        {
            iter ++;
            eqs->fitSP(i, tol, 1, 1, peak0);
            if (iter%5==0) eqs->fitSP(i, tol, 5, 0, peak0);
            eqs->calcSP(i);
            spi->fitError(Emax, Emean, CV);

            gr_post->setData(spi->qGetP(), spi->qGetSPA());
            qplot->rescaleAxes();
            qplot->xAxis->setRangeLower(0.01);
            qplot->yAxis->setRangeLower(0.0);
            qplot->replot();

            ui->progressBar->setValue(iter);
        }

        if (ui->progressBar->value()<mit)
            ui->progressBar->setValue(mit);
        ui->progressBar->hide();

        if (Emax <= 1.2*tol) {
            msg = tr("Spectrum Fitting Converged After %1 iterations!").arg(iter)
                    + tr("\nTotal Consumed Time: %1 s").arg((double)(clock()-ct)/CLOCKS_PER_SEC);
            QMessageBox::information(0, tr("EQSignal"), msg);
        }
        else {
            msg = tr("Spectrum Fitting not Converged After %1 iterations!").arg(iter)
                    + tr("\nTotal Consumed Time: %1 s").arg((double)(clock()-ct)/CLOCKS_PER_SEC);
            QMessageBox::warning(0, tr("EQSignal"), msg);
        }

    }

    ErrInfo = tr("Before Fitting")
            + tr("Max Error: %1%, Mean Error: %2%")
            .arg(Emaxp*100.0).arg(Emeanp*100.0)
            + ", " + tr("CV: %1").arg(CVp)
            + QString(";  ")
            + tr("After Fitting")
            + tr("Max Error: %1%, Mean Error: %2%")
            .arg(Emax*100.0).arg(Emean*100.0)
            + ", " + tr("CV: %1").arg(CV);
    ui->statusBar->showMessage(ErrInfo);

    eqs->a2vd();
    eqs0->copyAccFrom(eqs);
    eqs0->a2vd();
    this->plotTH(false);

}

void MainWindow::on_Paras_currentChanged(int index)
{
    switch (index) {
    case 7:
        ui->ParaTable->setItem(0,1,new QTableWidgetItem(QString::number(eqs->getPeakAcc())));
        ui->ParaTable->setItem(1,1,new QTableWidgetItem(QString::number(eqs->getPeakVel())));
        ui->ParaTable->setItem(2,1,new QTableWidgetItem(QString::number(eqs->getPeakDisp())));
        ui->ParaTable->setItem(3,1,new QTableWidgetItem(QString::number(eqs->getRMSAcc())));
        ui->ParaTable->setItem(4,1,new QTableWidgetItem(QString::number(eqs->getRMSVel())));
        ui->ParaTable->setItem(5,1,new QTableWidgetItem(QString::number(eqs->getRMSDisp())));
        break;
    default:
        break;
    }
}

void MainWindow::on_SPComPare_clicked()
{
    int i = ui->CDR->currentIndex();
    eqs->calcSP(i,false);

    ui->ViewSPA->clearGraphs();
    ui->tabWidget->setCurrentIndex(2);
    this->plotSPAi();

    double Emax, Emean, CV;
    eqs->getSP(i)->fitError(Emax, Emean, CV);

    FIT_SPA = true;

    QString ErrInfo = tr("Before Fitting")
            + tr("Max Error: %1%, Mean Error: %2%")
            .arg(Emax*100.0).arg(Emean*100.0)
            + ", " + tr("CV: %1").arg(CV);
    ui->statusBar->showMessage(ErrInfo);
}

void MainWindow::on_CalcRes_clicked()
{
    double zeta = ui->RDR->value();
    double P = ui->RPD->value();
    int method;

    double *cp = zeros(8);

    QString extra_para;
    QStringList extra_para_list;

    if (ui->NLR->isChecked()) {
        method = ui->NLMethod->currentIndex();
		cp[0] = ui->MU->value();
        cp[7] = (double)(ui->NLModel->currentIndex());

        extra_para = ui->ExtraPara->text().trimmed();
        extra_para_list = extra_para.split(",",QString::SkipEmptyParts);

        if (extra_para_list.count()>0) {
            cp[1] = extra_para_list[0].toDouble();
        }
		if (extra_para_list.count()>1) {
			cp[2] = extra_para_list[1].toDouble();
		}

        eqs->responseNL(zeta,P,method,cp);
    }
    else {
        method = ui->RSM->currentIndex();
        eqs->response(zeta,P,method);
    }

    plotRES();
    plotEnergy();
    plotHyst();
}

void MainWindow::on_actionGenWave_triggered()
{
//    gwdialog->show();
}

void MainWindow::on_Resample_clicked()
{
    eqs->resample(ui->RSR->value());
    eqs0->resample(ui->RSR->value());

    ui->TS->setValue(eqs->getDt()*2.5);

    plotTH();
}

void MainWindow::on_actionBasicInfo_triggered()
{
    QString msg;
    msg = tr("Name: %1\nTime Interval: %2\nNo. of Points: %3")
            .arg(eqsName).arg(eqs->getDt()).arg(eqs->getN());

    QMessageBox::information(0,"EQSignal",msg);
}

void MainWindow::on_tabWidget_currentChanged(int index)
{
    QString msg;
    double dt;
    int IPA, IPV, IPD;
    switch (index) {
    case 0:
        IPA = peakLoc(eqs->getAcc(),eqs->getN());
        IPV = peakLoc(eqs->getVel(),eqs->getN());
        IPD = peakLoc(eqs->getDisp(),eqs->getN());
        dt = eqs->getDt();

        msg = QString("PA=%1@%4  PV=%2@%5  PD=%3@%6")
                .arg(eqs->getAcc()[IPA])
                .arg(eqs->getVel()[IPV])
                .arg(eqs->getDisp()[IPD])
                .arg(dt*(IPA)).arg(dt*(IPV)).arg(dt*(IPD));
        ui->statusBar->showMessage(msg);
        break;
    case 2:
        if (FIT_SPA) {
            showSPAErrorMsg();
        }
        else
            ui->statusBar->clearMessage();
        break;
    default:
        ui->statusBar->clearMessage();
        break;
    }
}

void MainWindow::on_actionCalcSPA_triggered()
{
    setupSP();
    eqs->calcSP(true);
    plotSP();
    plotSPA();
    ui->Paras->setCurrentWidget(ui->SPA);
    FIT_SPA = false;
}

void MainWindow::on_actionFFT_triggered()
{
    on_CalcFFT_clicked();
    ui->Paras->setCurrentWidget(ui->FFT);
}

void MainWindow::on_actionFitSPA_triggered()
{
    ui->Paras->setCurrentWidget(ui->SPF);
    double Tg = QInputDialog::getDouble(this, tr("Characteristic Period"), "Tg = ",
                                        0.9, 0.0, 1.1, 1);
    double zeta = QInputDialog::getDouble(this, tr("Damping Ratio"), "Zeta = ",
                                        0.05, 0.0, 1.0, 2);
    ui->Tg->setValue(Tg);
    setupSP();

    int zi = ui->CDR->findText(QString::number(zeta));
    int zc = ui->CDR->count();
    if (zi>-1)
        ui->CDR->setCurrentIndex(zi);
    else {
        ui->CDR->addItem(QString::number(zeta));
        ui->CDR->setCurrentIndex(zc);
    }


    double PAF = 2.25;
    double SF = 1.0;

    eqs->setSPT(Tg,PAF,SF);
    eqs->calcSP();

    QStringList fms;
    fms << ui->FitMethod->itemText(0) << ui->FitMethod->itemText(1) << ui->FitMethod->itemText(2);

    QString fm = QInputDialog::getItem(this,"EQSignal",tr("FM"),fms,2,false);

    ui->FitMethod->setCurrentText(fm);

    QMessageBox::information(this,"EQSignal",tr("Press OK to Perform Spectrum Fitting."));

    on_SPFit_clicked();

}

void MainWindow::on_actionOpenLast_triggered()
{
    QString fileName = lastFile;

    if (fileName.isEmpty()) return;

    bool isTxt = fileName.endsWith(".txt") || (fileName.endsWith(".TXT"));
    bool isAt2 = fileName.endsWith(".at2") || (fileName.endsWith(".AT2"));

    QFileInfo fi(fileName);
    eqsName = fi.baseName();
    this->setWindowTitle("EQSignal -- " + eqsName);

    if (isTxt)
    {
        ui->URL->setText(fileName);
        double DT = QInputDialog::getDouble(this, eqsName + "'s " + tr("Time Interval"), "dt = ", 0.02, 0.0, 10.0, 3);
        this->readtxt(fileName, DT);

    }
    else if (isAt2)
    {
        ui->URL->setText(fileName);
        this->readnga(fileName);
    }
}

void MainWindow::on_actionEndtoZero_triggered()
{
    int N = eqs->getN();
    double pka = peak(eqs->getAcc(),N);
    double pkv = peak(eqs->getVel(),N);
    double pkd = peak(eqs->getDisp(),N);
    double aend = fabs(eqs->getAcc()[N-1]/pka);
    double vend = fabs(eqs->getVel()[N-1]/pkv);
    double dend = fabs(eqs->getDisp()[N-1]/pkd);
    if ( aend < 1.0e-3 && vend < 1.0e-3 && dend < 1.0e-3) return;

    int ntp = ui->NTP->value();
    eqs->endAlign(ntp,true,2);
    plotTH();
}
