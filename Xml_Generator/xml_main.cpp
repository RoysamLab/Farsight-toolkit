#include <QtGui/QApplication>
#include <QtCore/QObject>

 #include "Xml_Generator.h"

 int main(int argc, char *argv[])
 {
	QApplication app(argc, argv);


    QCoreApplication::setOrganizationName("RPI");
    QCoreApplication::setOrganizationDomain("farsight-toolkit.org");
    QCoreApplication::setApplicationName("Farsight Nucleus Editor");

     Xml_Generator xmlGen;
     xmlGen.show();
     return app.exec();
 }