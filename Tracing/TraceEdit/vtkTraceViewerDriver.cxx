/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#include <iostream>

#include "TraceView3D.h"
#include <QtGui/QApplication>
#include <QtCore/QObject>

#include "vtkCellArray.h"
#include "vtkPolyData.h"

int main (int argc, char* argv[])
{
	//argc will be taken care of in QT app
	QApplication app(argc, argv);
	app.setOrganizationName("FARSIGHT Toolkit");
	app.setOrganizationDomain("www.farsight-toolkit.org");
	app.setApplicationName("Trace Editor");
	app.setApplicationVersion("V2.0");
	View3D *View = new View3D();
	View->show();
  int retval = View->runTests();
  if(retval == -1)
  {
    retval = app.exec();
  }
  else
  {
    View->close();
  }
	delete View;
  std::cout << "TraceEdit returns " << retval << std::endl;
	return retval;
}
