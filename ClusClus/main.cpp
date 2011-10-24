#include <QtGui/QApplication>
#include "ClusClusMainwindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    ClusClusMainWindow ccw;
    ccw.show();

    return app.exec();
}