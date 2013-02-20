/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

//This is the main GUI application.  It consists of a control bar and any number of windows
// similar to the ImageJ setup.
//The controlBar is the parent widget of all of the other windows and contains all of the shared models
//The windows provide various views into underlying data
//Simple image viewing capability is also possible.

#include <QtGui/QApplication>
#include <QtCore/QObject>

#include "NucleusEditor.h"

int main(int argc, char *argv[])
{
    	QApplication app(argc, argv);

	QCoreApplication::setOrganizationName("RPI");
    	QCoreApplication::setOrganizationDomain("farsight-toolkit.org");
    	QCoreApplication::setApplicationName("Farsight Nucleus Editor");

	NucleusEditor nucEdit(0,0);
	nucEdit.show();

  int retval = nucEdit.runTest();
  if(retval == -1)
  {
    retval = app.exec();
  }
  return retval;
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
