#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <QInputDialog>
#include <QTextBrowser>
#include <QTableWidget>
#include <QComboBox>
#include <QPen>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QSpinBox>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QList>
#include <QString>
#include <QMimeData>
#include <QDialogButtonBox>
#include <QDialog>
#include "EQSignal.h"
#include "eqtablewidget.h"
#include "qcustomplot.h"

namespace Ui {
class MainWindow;
}

class GenWaveDialog : public QDialog
{
	Q_OBJECT
public:
	GenWaveDialog(QWidget *parent = 0) : QDialog(parent)
	{
		this->resize(400, 300);
		gridLayout = new QGridLayout(this);

		label_wavetype = new QLabel(tr("Wave Type"), this);
		label_N = new QLabel("N", this);
		label_dt = new QLabel("dt", this);
		label_A0 = new QLabel("A0", this);
		label_A = new QLabel("A", this);
		label_T = new QLabel("T", this);
		label_phi = new QLabel("phi", this);

		wavetype = new QComboBox(this);
        wavetype->addItem(tr("A0 + A*sin(2*pi/T*t + phi)"));
		N = new QSpinBox(this);
        N->setRange(128,65536);
        N->setValue(2048);
		dt = new QDoubleSpinBox(this);
        dt->setRange(0.001,0.1);
        dt->setValue(0.02);
		A0 = new QDoubleSpinBox(this);
		A = new QDoubleSpinBox(this);
        A->setValue(1.0);
		T = new QDoubleSpinBox(this);
        T->setValue(1.0);
		phi = new QDoubleSpinBox(this);

		gridLayout->addWidget(label_wavetype, 0, 0, 1, 1);
		gridLayout->addWidget(label_N, 1, 0, 1, 1);
		gridLayout->addWidget(label_dt, 2, 0, 1, 1);
		gridLayout->addWidget(label_A0, 3, 0, 1, 1);
		gridLayout->addWidget(label_A, 1, 2, 1, 1);
		gridLayout->addWidget(label_T, 2, 2, 1, 1);
		gridLayout->addWidget(label_phi, 3, 2, 1, 1);
		gridLayout->addWidget(wavetype, 0, 1, 1, 3);
		gridLayout->addWidget(N, 1, 1, 1, 1);
		gridLayout->addWidget(dt, 2, 1, 1, 1);
		gridLayout->addWidget(A0, 3, 1, 1, 1);
		gridLayout->addWidget(A, 1, 3, 1, 1);
		gridLayout->addWidget(T, 2, 3, 1, 1);
		gridLayout->addWidget(phi, 3, 3, 1, 1);

		buttonBox = new QDialogButtonBox(this);
		buttonBox->setOrientation(Qt::Horizontal);
		buttonBox->setStandardButtons(QDialogButtonBox::Cancel | QDialogButtonBox::Ok);

		gridLayout->addWidget(buttonBox, 4, 0, 1, 4);

		connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
		connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
	}

	QLabel *label_wavetype;
	QLabel *label_N;
	QLabel *label_dt;
	QLabel *label_A0;
	QLabel *label_A;
	QLabel *label_T;
	QLabel *label_phi;
	QComboBox *wavetype;
	QDoubleSpinBox *dt;
	QDoubleSpinBox *A0;
	QDoubleSpinBox *A;
    QDoubleSpinBox *T;
	QDoubleSpinBox *phi;
	QSpinBox *N;

	QGridLayout *gridLayout;
	QDialogButtonBox *buttonBox;

};

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

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void readtxt(const char *filename, double dt);
    void readtxt(QString filename, double dt);
    void readnga(const char *filename);
    void readnga(QString filename);

    static QPen autoPen(int i);
    static QColor reverseColor(QColor c);

protected:
    void dragEnterEvent(QDragEnterEvent *event) {
        if (event->mimeData()->hasFormat("text/uri-list"))
            event->acceptProposedAction();
    }

    void dropEvent(QDropEvent *event);

public slots:

    void showDriftMsg();
    void confirm();
    void saveTH();
    void saveSP();
    void setupSP();

    void plotTH(bool changeTab=true);
    void plotRES();
    void plotEnergy();
    void plotHyst();
    void plotSPA();
    void plotSP();
    void plotFFT();
    void plotPSD();
    void plotSPT();
    void plotSPAi();

    void SPXScaleChanged();

    void fillSPTable(int k=0);
    void fillTHTable();

    void pdw_applyButton_clicked();
    void sptw_applyButton_clicked();

    void fillRESTable();
	void genWave();

private slots:

    void on_Norm_clicked();

    void on_Reload_clicked();

    void on_Filt_clicked();

    void on_Detrend_clicked();

    void on_Align_clicked();

    void on_Adjust_clicked();

    void on_SetTrimEdge_clicked();

    void on_Trim_clicked();

    void on_CalcFFT_clicked();

    void on_CalcPSD_clicked();

    void on_CalcSPA_clicked();

    void on_CalcSP_clicked();

    void on_DR_textChanged(const QString &arg1);

    void on_actionOpen_triggered();

    void on_actionSaveData_triggered();

    void on_actionSaveFig_triggered();

    void on_actionSaveAcc_triggered();

    void on_GenSPT_clicked();

    void on_actionViewData_triggered();

    void on_DM_activated(int index);

    void on_TSPT_activated(int index);

    void on_SPFit_clicked();
    
    void on_Paras_currentChanged(int index);

    void on_SPComPare_clicked();

    void on_CalcRes_clicked();

    void on_actionGenWave_triggered();

    void on_Resample_clicked();

    void on_actionBasicInfo_triggered();

    void on_tabWidget_currentChanged(int index);

    void on_actionCalcSPA_triggered();

    void on_actionFFT_triggered();

    void on_actionFitSPA_triggered();

private:
    void setupConnections();
	void initViewTH();
    void initViewRES();
    void initViewHyst();
	void initViewSPA();
    void initViewSP();
    void initViewFFT();
    void initViewEnergy();
    void initTable();

    void clearView(QCustomPlot *qplot);

    void readConfig();
    void writeConfig();

    Ui::MainWindow *ui;

    GenWaveDialog *gwd;

    EQSignal *eqs;
    EQSignal *eqs0;

    QDir workDir;
    QString eqsName;

    bool saveAccWithTime;
    bool normOnRead;

    EQTableWidget *dataTableTH;
    EQTableWidget *dataTableRES;
    EQTableWidget *dataTableSP;

    PeriodsDefineWidget *pdw;
    SPTDefineWidget *sptw;

    QWidget *wth;
    QWidget *wres;
    QWidget *wsp;

    QComboBox *dr;

    void ViewTHData();
    void ViewSPData();

};

#endif // MAINWINDOW_H
