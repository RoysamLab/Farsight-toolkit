#include <iostream>

#include "CurveletGUI.h"
#include <QtGui/QApplication>
#include <QtCore/QObject>
int main (int argc, char* argv[])
{
	//argc will be taken care of in QT app
	QApplication app(argc, argv);
	app.setOrganizationName("FARSIGHT Toolkit");
	app.setOrganizationDomain("www.farsight-toolkit.org");
	app.setApplicationName("Curvelet Gui");
	app.setApplicationVersion("V1.0");
	CurveletGUI * mainCurveletGUI=new CurveletGUI();
	mainCurveletGUI->show();
	int retval = app.exec();
	delete mainCurveletGUI;
	return retval;
}
