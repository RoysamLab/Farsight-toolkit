Initial cache files for CMake configuration of various dependencies.

1) Install prerequisites (Boost, QT, GLUT) for your platform.
  * On Debian/Ubuntu, use packages: libboost-dev, qt4-qmake, libglut3-dev
  * On Mac OS X, use MacPorts packages: boost, qt4-mac, mesa

For more details, see:

  http://farsight-toolkit.org/wiki/FARSIGHT_HowToBuild

2) Download source code for VXL, VTK and ITK (either tarballs or from CVS).

3) Edit the configuration file according to your directory structure:

  vi vars.conf

4) Compile everything with:

  bash build-all.sh


============================
Known working configurations
============================

   OS: Mac OS X 10.5.8 with Xcode 3.1.3
Boost: Installed boost package using MacPorts on 2009-10-24.
   QT: Installed qt4-mac package using MacPorts on 2009-10-24.
 GLUT: Installed mesa package using MacPorts on 2009-10-24.
  VXL: 1.13.0 release version, compiled from source
  ITK: 3.16.0 release version, compiled from source
  VTK: 5.4.2 release version, compiled from source
  FTK: r1146 from Subversion
