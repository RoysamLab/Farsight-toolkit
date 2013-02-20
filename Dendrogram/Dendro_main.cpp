#include <QtGui/QApplication>
#include <QtCore/QObject>
#include "Dendrogram.h"
 int main(int argc, char *argv[])
 {
	QApplication app(argc, argv);

	
	
	
 
	
	/*QCoreApplication::setOrganizationName("RPI");
    QCoreApplication::setOrganizationDomain("farsight-toolkit.org");
    QCoreApplication::setApplicationName("Farsight Nucleus Editor");*/

     Dendrogram Dendro_Selection(0,0);
	 //Dendro_Selection.show();
	
	
	return app.exec();
	// return 1;
 }
