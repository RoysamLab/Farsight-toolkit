This folder of CMake configuration cache files and build scripts makes
it much simpler to build FARSIGHT on Mac OS X and Linux systems.

Steps:

1) Install prerequisites (Boost, QT, GLUT) for your platform.
  * On Mac OS X, use MacPorts packages: boost, qt4-mac, mesa
  * On Debian/Ubuntu, use packages: libboost-dev, qt4-qmake, libglut3-dev

For more details, see:

  http://farsight-toolkit.org/wiki/FARSIGHT_HowToBuild

2) Download source code for VXL, VTK and ITK
   (either tarballs or from version control).

3) Edit the configuration file according to your directory structure:

  vi vars.conf

4) Compile everything with:

  bash build-all.sh


============================
Known working configurations
============================

   OS: Mac OS X 10.6.4
Boost: Installed boost package using MacPorts on 2010-08-17
   QT: Installed qt4-mac package using MacPorts on 2010-08-17
 GLUT: Installed mesa package using MacPorts on 2010-08-17
  VXL: 1.13.0 release version, compiled from source
  ITK: 3.20.0 release version, compiled from source
  VTK: 5.6.0 release version, compiled from source
  FTK: r2066 from Subversion
