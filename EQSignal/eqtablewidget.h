#ifndef EQTABLEWIDGET_H
#define EQTABLEWIDGET_H

#include <QTableWidget>


class EQTableWidget : public QTableWidget
{
    Q_OBJECT
public:
    EQTableWidget(QWidget *parent=0) : QTableWidget(parent) {}

    double *getColumnData(int col);
    double *getRowData(int row);

    void setColumn(int col, double *data);
    void setColumn(int col, QStringList data);

signals:
    void rowNumberChanged(int);

public slots:
    void setColumnNumber(int nc) {this->setColumnCount(nc);}
    void setRowNumber(int nr) {this->setRowCount(nr);}

protected:
    void keyPressEvent(QKeyEvent *e);

};

#endif // EQTABLEWIDGET_H
