#include <iostream>

#include "Seed3D.h"
#include "Seed3D.cxx"
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
	Seed3D Seeds(argc, argv);
  Seeds.show();
  return app.exec();
  }
