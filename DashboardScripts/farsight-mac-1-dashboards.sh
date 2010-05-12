#!/bin/sh
export DISPLAY=:0.0
/usr/local/bin/ctest -V -VV -S /projects/Dashboards/farsight-trunk/DashboardScripts/farsight-mac-1_nightly_release_makefiles.cmake > /projects/Dashboards/Logs/nightly_release_makefiles.log 2>&1
/usr/local/bin/ctest -V -VV -S /projects/Dashboards/farsight-trunk/DashboardScripts/farsight-mac-1_nightly_debug_makefiles.cmake > /projects/Dashboards/Logs/nightly_debug_makefiles.log 2>&1
/usr/local/bin/ctest -V -VV -S /projects/Dashboards/farsight-trunk/DashboardScripts/farsight-mac-1_continuous_debug_makefiles.cmake > /projects/Dashboards/Logs/continuous_debug_makefiles.log 2>&1

