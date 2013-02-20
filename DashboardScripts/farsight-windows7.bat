call "C:\Program Files (x86)\CMake 2.8\bin\ctest.exe" -S C:/dashboard/src/farsight/DashboardScripts/farsight-windows7_dbg_static_nightly.cmake -VV -O C:/dashboard/farsight_dbg_static_nightly_log.txt
cd C:/dashboard/src/boost
svn up
call "C:\Program Files (x86)\CMake 2.8\bin\ctest.exe" -S C:/dashboard/src/farsight/DashboardScripts/farsight-windows7_vxl_dbg_nightly.cmake -VV -O C:/dashboard/vxl_dbg_nightly_log.txt
call "C:\Program Files (x86)\CMake 2.8\bin\ctest.exe" -S C:/dashboard/src/farsight/DashboardScripts/farsight-windows7_itk_dbg_nightly.cmake -VV -O C:/dashboard/itk_dbg_nightly_log.txt
call "C:\Program Files (x86)\CMake 2.8\bin\ctest.exe" -S C:/dashboard/src/farsight/DashboardScripts/farsight-windows7_vtk_dbg_nightly.cmake -VV -O C:/dashboard/vtk_dbg_nightly_log.txt
call "C:\Program Files (x86)\CMake 2.8\bin\ctest.exe" -S C:/dashboard/src/farsight/DashboardScripts/farsight-windows7_dbg_nightly.cmake -VV -O C:/dashboard/farsight_dbg_nightly_log.txt