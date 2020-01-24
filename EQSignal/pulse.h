#ifndef PULSE_H
#define PULSE_H

#include <QDialog>

namespace Ui {
class Pulse;
}

class Pulse : public QDialog
{
    Q_OBJECT

public:
    explicit Pulse(QWidget *parent = nullptr);
    ~Pulse();

	void plot(double *t, double *v, int n);
    int type();
	double Tp();
	double t0();
	double t1();
	double gamma();
	double vp();

	QPushButton *genButton();
	QPushButton *addButton();

private:
    Ui::Pulse *ui;
};

#endif // PULSE_H
