#include "pulse.h"
#include "ui_pulse.h"
#include "eqs.h"

Pulse::Pulse(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Pulse)
{
    ui->setupUi(this);
	ui->pulseView->clearGraphs();
}

Pulse::~Pulse()
{
    delete ui;
}

void Pulse::plot(double * t, double * v, int n)
{
	QVector<double> Time, Vel;
	Time = A2QV(t, n);
	Vel = A2QV(v, n);

	QCustomPlot *qplot = ui->pulseView;
	qplot->clearGraphs();
	QCPGraph *gr = qplot->addGraph();
	gr->setData(Time, Vel);
	qplot->rescaleAxes();

	qplot->xAxis->setLabel(tr("Time"));
	qplot->yAxis->setLabel(tr("Velocity"));

	qplot->replot();
}

int Pulse::type()
{
	return ui->type->currentIndex();
}

double Pulse::Tp()
{
	return ui->Tp->value();
}

double Pulse::t0()
{
	return ui->t0->value();
}

double Pulse::t1()
{
	return ui->t1->value();
}

double Pulse::gamma()
{
	return ui->gamma->value();
}

double Pulse::vp()
{
	return ui->vp->value();
}

QPushButton * Pulse::genButton()
{
	return ui->genPulse;
}

QPushButton * Pulse::addButton()
{
	return ui->addPulse;
}
