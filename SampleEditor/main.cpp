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

//This is the main GUI application.  It consists of a control bar and any number of windows
// similar to the ImageJ setup.
//The controlBar is the parent widget of all of the other windows and contains all of the shared models
//The windows provide various views into underlying data
//Simple image viewing capability is also possible.

#include <QtGui/QApplication>
#include <QtCore/QObject>

#include "SampleEditor.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

	SampleEditor sampleEdit(0,0);
	sampleEdit.show();

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
