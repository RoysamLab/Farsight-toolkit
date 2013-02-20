#include <QtGui/QApplication>
#include "clusclusMainwindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    ClusClusMainWindow ccw;
    ccw.show();

    return app.exec();
}
