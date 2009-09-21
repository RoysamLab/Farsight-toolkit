#include <iostream>

#include "Seed3D.h"

#include <QtGui/QApplication>
#include <QtCore/QObject>

// to bypass linking error 2019 
// 2019 unresolved external symbol _WinMain@16
// go to project properties->Linker -> System -> SubSystem: choose Console
int main (int argc, char  **argv)
  {
 
  QApplication app(argc, argv);
 
	Seed3D Seeds(0, 0);
	
  Seeds.show();
  
  return app.exec();
  }

//This function is used to create a GUI application that does not show the console.
// See SET_TARGET_PROPERTIES(Farsight PROPERTIES WIN32_EXECUTABLE 0)
// in CMakeLists.txt to change the property.
#ifdef _WIN32
#include <windows.h>
int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
		     LPSTR lpCmdLine, int nCmdShow)
{
  return main(nCmdShow, &lpCmdLine);
}
#endif
