#include <iostream>

#include "View3D.h"
#include <QtGui/QApplication>
#include <QtCore/QObject>

int main (int argc, char* argv[])
  {
  if(argc < 2)
    {
    cerr << argv[0] << " <filename>" << endl;
    return 0;
    }
  QApplication app(argc, argv);
	View3D View(argc, argv);
  View.show();
  return app.exec();
  }
