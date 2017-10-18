#ifndef GENWAVEDIALOG
#define GENWAVEDIALOG

#include <QtCore>
#include <QDialog>
#include <QGridLayout>
#include <QLabel>
#include <QComboBox>
#include <QSpinBox>
#include <QDialogButtonBox>

class GenWaveDialog : public QDialog
{
    Q_OBJECT
public:
    GenWaveDialog(QWidget *parent = 0);

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

private slots:
	void on_wavetype_activated(int index);
};

#endif // GENWAVEDIALOG

