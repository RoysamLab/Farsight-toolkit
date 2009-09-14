Initial cache files for CMake configuration of various dependencies.

For example, to compile VXL as needed for FARSIGHT, use:

  cmake -C /path/to/farsight/cache/VXL.cmake /path/to/VXL

Or edit the configuration file and then use the shell scripts:

  vi vars.conf
  bash build-vxl.sh
  bash build-itk.sh
  bash build-vtk.sh
  bash build-ftk.sh

Or compile everything:

  vi vars.conf
  bash build-all.sh


=========
VTK notes
=========

After installing VTK, even with VTK_USE_GUISUPPORT and VTK_USE_QVTK set to ON,
two VTK header files will be missing from the installation directory
(/usr/local/include/vtk-5.4):

  GUISupport/Qt/vtkQtBarChartView.h
  GUISupport/Qt/vtkQtChartViewBase.h

The build-vtk.sh script will try to copy these files for you.


=========
FTK notes
=========

To install GLUT on Linux, you can use freeglut:

  sudo aptitude install libglut3-dev
