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

#include "View3D.h"
#include <QtGui/QApplication>
#include <QtCore/QObject>

#include "vtkCellArray.h"
#include "vtkPolyData.h"

int main (int argc, char* argv[])
  {
  //if(argc < 2)
  //  {
  //  cerr << argv[0] << " <filename>" << endl;
  //  return 0;
  //  }
  QApplication app(argc, argv);
  app.setOrganizationName("FARSIGHT Toolkit");
  app.setOrganizationDomain("farsight-toolkit.org");
  app.setApplicationName("Trace Editor");
	View3D *View = new View3D(argc, argv);
  View->show();
  int retval = app.exec();
  delete View;
  return retval;
  }
