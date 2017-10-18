#include "mainwindow.h"
#include <QApplication>
#include <QSplashScreen>
#include <QDesktopWidget>
#include <QTranslator>
#include <QString>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QString qm("/trans.qm");
    qm.prepend(QApplication::applicationDirPath());
    QTranslator *trans = new QTranslator;
    trans->load(qm);
    a.installTranslator(trans);

    QSplashScreen *splash = new QSplashScreen;
    splash->setPixmap(QPixmap(":/MainWindow/splash.png"));
    splash->show();

    QFont font;
    font.setFamily("Microsoft YaHei");
    a.setFont(font);
    MainWindow w;

    w.move(320, 100);

    w.show();
    splash->finish(&w);

    return a.exec();
}
