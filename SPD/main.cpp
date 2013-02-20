#include <QtGui/QApplication>
#include "spdmainwindow.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    SPDMainWindow w;
    w.show();

    return a.exec();
}
