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
#include "spectradefinewidget.h"
#include "genwavedialog.h"
#include "qcustomplot.h"

namespace Ui {
class MainWindow;
}

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

    void closeEvent(QCloseEvent *event)
    {
        writeConfig();
        event->accept();
    }

signals:
    void xsChanged();

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
    void plotSPABeforeFit();
    void plotSPAAfterFit();

    void setXScale();

    void fillSPTable(int k=0);
    void fillTHTable();

    void pdw_applyButton_clicked();
    void sptw_applyButton_clicked();

    void fillRESTable();
	void genWave();
    void showSPAErrorMsg();

    void genArtificialEQWave(double *a, int N, double DT);
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

//    void on_DR_textChanged(const QString &arg1);

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

    void on_actionOpenLast_triggered();

    void on_actionEndtoZero_triggered();

    void on_actionChangSpecialPoint_triggered();

    void on_AccEndsToZero_clicked();

    void on_actionToggleXScale_triggered();

    void on_Recover_clicked();

    void on_DSPA_clicked();

    void on_DFit_clicked();

    void on_GenAW_clicked();

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
    QString lastFile;
    QString eqsName;

    bool saveAccWithTime;
    bool normOnRead;
    bool XS_LOG;

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

    bool FIT_SPA; // Control whether to show spectrum fitting information

};

#endif // MAINWINDOW_H
