#!/bin/sh
export DISPLAY=:0.0
/usr/local/bin/ctest -V -VV -S /Dashboards/farsight-trunk-nightly/DashboardScripts/farsight-ubuntu-1_gcc43_rel_static_nightly.cmake > /Dashboards/Logs/rel_static_nightly.log 2>&1
/usr/local/bin/ctest -V -VV -S /Dashboards/farsight-trunk-nightly/DashboardScripts/farsight-ubuntu-1_gcc43_dbg_static_nightly.cmake > /Dashboards/Logs/dbg_static_nightly.log 2>&1
/usr/local/bin/ctest -V -VV -S /Dashboards/farsight-trunk-nightly/DashboardScripts/farsight-ubuntu-1_gcc43_dbg_static_nightly_vtkgit.cmake > /Dashboards/Logs/dbg_static_nightly.log 2>&1
/usr/local/bin/ctest -V -VV -S /Dashboards/farsight-trunk-continuous/DashboardScripts/farsight-ubuntu-1_gcc43_dbg_static_continuous.cmake > /Dashboards/Logs/dbg_static_continuous.log 2>&1
