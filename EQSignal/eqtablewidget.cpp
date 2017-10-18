#include "eqtablewidget.h"
#include <QApplication>
#include <QKeyEvent>
#include <QClipboard>

double *EQTableWidget::getColumnData(int col)
{
    int nr = this->rowCount();
    double *data = new double[nr];

    QTableWidgetItem *item;

    for (int i=0; i<nr; ++i) {
        item = this->item(i,col);
        if (0 == item)
            data[i] = 0.0;
        else
            data[i] = item->text().toDouble();
    }
    return data;
}

double *EQTableWidget::getRowData(int row)
{
    int nc = this->columnCount();
    double *data = new double[nc];

    QTableWidgetItem *item;

    for (int i=0; i<nc; ++i) {
        item = this->item(row,i);
        if (0 == item)
            data[i] = 0.0;
        else
            data[i] = item->text().toDouble();
    }
    return data;
}

void EQTableWidget::setColumn(int col, double *data)
{
    int nr = rowCount();

    for (int i=0; i<nr; ++i) {
		setItem(i, col, new QTableWidgetItem(QString::number(data[i])));
    }
}

void EQTableWidget::setColumn(int col, QStringList data)
{
    int nr = rowCount();

    for (int i=0; i<nr; ++i) {
        setItem(i, col, new QTableWidgetItem(data[i]));
    }
}

void EQTableWidget::keyPressEvent(QKeyEvent *e)
{
    int rt,rb,cl,cr;
    QTableWidgetSelectionRange range;
    QString s, line;
    QStringList lines, cols;
    QClipboard *clipboard = QApplication::clipboard();

    if (e->modifiers() & Qt::ControlModifier) {
        switch (e->key()) {
        case Qt::Key_C:
            range = this->selectedRanges()[0];
            rt = range.topRow();
            rb = range.bottomRow();
            cl = range.leftColumn();
            cr = range.rightColumn();

            for (int i=rt; i<=rb; ++i) {
                for (int j=cl; j<cr; ++j) {
                    s.append(this->item(i,j)->text() + "\t");
                }
                s.append(this->item(i,cr)->text() + "\n");
            }

            clipboard->setText(s);
            break;
        case Qt::Key_V:
            range = this->selectedRanges()[0];
            rt = range.topRow();
            rb = range.bottomRow();
            cl = range.leftColumn();
            cr = range.rightColumn();

            s = clipboard->text();
            s = s.remove("\r");

            lines = s.split("\n",QString::SkipEmptyParts);

            line = lines[0];
            cols = line.split("\t",QString::SkipEmptyParts);

            if ( (rt + lines.count()) > this->rowCount() ) {
                this->setRowCount(rt + lines.count());
                emit rowNumberChanged(this->rowCount());
            }
            if ( (cl + cols.count()) > this->columnCount() ) this->setColumnCount(cl + cols.count());

            for (int i=rt; i<(rt + lines.count()); ++i) {
                line = lines[i-rt];
                cols = line.split("\t",QString::SkipEmptyParts);
                for (int j=cl; j<(cl + cols.count()); ++j) {
                    this->setItem(i,j,new QTableWidgetItem(cols[j-cl]));
                }
            }
            break;
        case Qt::Key_A:
            this->selectAll();
            break;
        default:
            break;
        }
    }
    else {
        QTableWidget::keyPressEvent(e);
    }
}

