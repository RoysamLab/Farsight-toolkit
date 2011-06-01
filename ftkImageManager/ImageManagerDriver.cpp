#include <iostream>

#include "ImageManager.h"
#include <QtGui/QApplication>
#include <QtCore/QObject>

int main (int argc, char* argv[])
{
	//argc will be taken care of in QT app
	QApplication app(argc, argv);
	app.setOrganizationName("FARSIGHT Toolkit");
	app.setOrganizationDomain("www.farsight-toolkit.org");
	app.setApplicationName("Image Manager");
	app.setApplicationVersion("V1.0");
	ImageFileManger * manager = new ImageFileManger();
	manager->show();
	int retval = app.exec();
	delete manager;
	return retval;
}
