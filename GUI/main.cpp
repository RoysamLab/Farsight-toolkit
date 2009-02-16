//This is the main GUI application.  It consists of a control bar and any number of windows
// similar to the ImageJ setup.
//The controlBar is the parent widget of all of the other windows and contains all of the shared models
//The windows provide various views into underlying data
//Simple image viewing capability is also possible.

#include <QtGui/QApplication>
#include <QtCore/QObject>

#include "ControlBar.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
	//Q_INIT_RESOURCE(farsight);

	ControlBar cbar(argv[0]);
	cbar.show();

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
