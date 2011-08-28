#!/bin/sh
export DISPLAY=:0.0
cd /home/gramak/dashboard/boost/
svn up
/usr/local/bin/ctest -V -VV -S /home/gramak/dashboard/farsight/DashboardScripts/farsight-opensuse-vxl_dbg_nightly.cmake > /home/gramak/dashboard/Logs/vxl_nightly_dbg.log 2>&1
/usr/local/bin/ctest -V -VV -S /home/gramak/dashboard/farsight/DashboardScripts/farsight-opensuse-itk_dbg_nightly.cmake > /home/gramak/dashboard/Logs/itk_nightly_dbg.log 2>&1
/usr/local/bin/ctest -V -VV -S /home/gramak/dashboard/farsight/DashboardScripts/farsight-opensuse-vtk_dbg_nightly.cmake > /home/gramak/dashboard/Logs/vtk_nightly_dbg.log 2>&1
/usr/local/bin/ctest -V -VV -S /home/gramak/dashboard/farsight/DashboardScripts/farsight-opensuse-gcc45_dbg_nightly.cmake > /home/gramak/dashboard/Logs/farsight_nightly_dbg.log 2>&1
