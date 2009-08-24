Initial cache files for CMake configuration of various dependencies.

For example, to compile VXL as needed for FARSIGHT, use:

  cmake -C /path/to/farsight/cache/VXL.cmake /path/to/VXL

Or use the shell scripts:

  VXL_DIR=/path/to/VXL bash build-vxl.sh
  ITK_DIR=/path/to/ITK bash build-itk.sh
  VTK_DIR=/path/to/VTK bash build-vtk.sh
  FTK_DIR=/path/to/FTK bash build-ftk.sh

Or edit the master shell script and then run to compile everything:

  vi build-all.sh
  bash build-all.sh

=========
VTK notes
=========

After installing VTK on Mac OS X, even though VTK_USE_GUISUPPORT and
VTK_USE_QVTK were both set to ON, needed to copy two missing header files
into the installation directory (/usr/local/include/vtk-5.4):

  GUISupport/Qt/vtkQtBarChartView.h
  GUISupport/Qt/vtkQtChartViewBase.h
