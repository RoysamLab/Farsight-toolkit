include_directories(${TraceEdit_SOURCE_DIR})

add_executable(TraceEdit-StartStop TraceEditStartStop.cxx)
target_link_libraries(TraceEdit-StartStop Trace)
add_test(TraceEdit-StartStop ${Farsight_BINARY_DIR}/exe/TraceEdit-StartStop)
set_tests_properties(TraceEdit-StartStop PROPERTIES
  FAIL_REGULAR_EXPRESSION "vtkDebugLeaks has detected LEAKS!"
)

if(QtTestingFound)

  if(APPLE)
    set(TraceEdit_EXE ${Farsight_BINARY_DIR}/exe/TraceEdit.app/Contents/MacOS/TraceEdit) 
    set(AppName "Farsight")
  else()
    set(TraceEdit_EXE ${Farsight_BINARY_DIR}/exe/TraceEdit)
    set(AppName "TraceEdit")
  endif()

  set(TraceEdit_TEST_DIR ${TraceEdit_SOURCE_DIR}/Testing)
  
  add_test(
    TraceEdit-ClearSettings
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C5TracedPoints.swc
    ${TraceEdit_TEST_DIR}/test_clearsettings.xml
  )
  set_tests_properties(TraceEdit-ClearSettings PROPERTIES
    PASS_REGULAR_EXPRESSION "All QSettings have been reverted to their default values"
  )

  add_test(
    TraceEdit-AutoMerge
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C5TracedPoints.swc
    ${TraceEdit_TEST_DIR}/test_automerge.xml
    ${TraceEdit_TEST_DIR}/baseline_automerge.png
  )
  set_tests_properties(TraceEdit-AutoMerge PROPERTIES
    DEPENDS TraceEdit-ClearSettings
  )

  add_test(
    TraceEdit-SelectTree
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/cnic_macaque_pyramidal.swc
    ${TraceEdit_TEST_DIR}/test_select_tree.xml
    ${TraceEdit_TEST_DIR}/baseline_select_tree.png
  )
  set_tests_properties(TraceEdit-SelectTree PROPERTIES
    DEPENDS TraceEdit-AutoMerge
  )

  add_test(
    TraceEdit-RaycastTools
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    ${TraceEdit_TEST_DIR}/test_raycast_tools.xml
    ${TraceEdit_TEST_DIR}/baseline_raycast_tools.png
  )
  set_tests_properties(TraceEdit-RaycastTools PROPERTIES
    DEPENDS TraceEdit-SelectTree
  )
    
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_load_trace_image_soma.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_load_trace_image_soma.xml
    @ONLY)
  add_test(
    TraceEdit-LoadTraceImageSoma
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${FARSIGHT_DATA_ROOT}/TraceData/
    ${CMAKE_CURRENT_BINARY_DIR}/test_load_trace_image_soma.xml
    ${TraceEdit_TEST_DIR}/baseline_load_trace_image_soma.png
  )
  set_tests_properties(TraceEdit-LoadTraceImageSoma PROPERTIES
    DEPENDS TraceEdit-RaycastTools
  )
  
  add_test(
    TraceEdit-Reload
    ${TraceEdit_EXE}
    ${TraceEdit_TEST_DIR}/test_reload.xml
    ${TraceEdit_TEST_DIR}/baseline_reload.png
    reload
  )
  set_tests_properties(TraceEdit-Reload PROPERTIES
    DEPENDS TraceEdit-LoadTraceImageSoma
  )
  
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_saveproject.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_saveproject.xml
    @ONLY)
  add_test(
    TraceEdit-SaveProject
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${CMAKE_CURRENT_BINARY_DIR}/test_saveproject.xml
    ${TraceEdit_TEST_DIR}/baseline_saveproject.png
  )
  set_tests_properties(TraceEdit-SaveProject PROPERTIES
    DEPENDS TraceEdit-Reload
    PASS_REGULAR_EXPRESSION "project saved to disk"
  )

  add_test(
    TraceEdit-LoadProject
    ${TraceEdit_EXE}
    ${TraceEdit_TEST_DIR}/test_loadproject.xml
    ${Farsight_BINARY_DIR}/Testing/Temporary/baseline_loadproject.png
    ${Farsight_BINARY_DIR}/Testing/Temporary/project_test.xml
  )
  set_tests_properties(TraceEdit-LoadProject PROPERTIES
    DEPENDS TraceEdit-SaveProject
  )
  
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_tracing.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_tracing.xml
    @ONLY)
  add_test(
    TraceEdit-Tracing
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    ${CMAKE_CURRENT_BINARY_DIR}/test_tracing.xml
    ${TraceEdit_TEST_DIR}/baseline_tracing.png
  )
  set_tests_properties(TraceEdit-Tracing PROPERTIES
    DEPENDS TraceEdit-LoadProject
    PASS_REGULAR_EXPRESSION "Starting Amit's Tracer"
  )
  
  add_test(
    TraceEdit-Settings
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_settings.xml
    ${TraceEdit_TEST_DIR}/baseline_settings.png
  )
  set_tests_properties(TraceEdit-Settings PROPERTIES
    DEPENDS TraceEdit-Tracing
  )
  
  add_test(
    TraceEdit-Gridlines
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_gridlines.xml
    ${TraceEdit_TEST_DIR}/baseline_gridlines.png
  )
  set_tests_properties(TraceEdit-Gridlines PROPERTIES
    DEPENDS TraceEdit-Settings
  )
  
  add_test(
    TraceEdit-MergeSelected
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C5TracedPoints.swc
    ${TraceEdit_TEST_DIR}/test_mergeselected.xml
    ${TraceEdit_TEST_DIR}/baseline_mergeselected.png
  )
  set_tests_properties(TraceEdit-MergeSelected PROPERTIES
    DEPENDS TraceEdit-Gridlines
  )
  
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_exportcells.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_exportcells.xml
    @ONLY)
  add_test(
    TraceEdit-ExportCells
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/cnic_macaque_pyramidal.swc
    ${CMAKE_CURRENT_BINARY_DIR}/test_exportcells.xml
    ${TraceEdit_TEST_DIR}/baseline_select_tree.png
  )
  set_tests_properties(TraceEdit-ExportCells PROPERTIES
    DEPENDS TraceEdit-Gridlines
    PASS_REGULAR_EXPRESSION "file written"
  )
  
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_raycast_somas.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_raycast_somas.xml
    @ONLY)
  add_test(
    TraceEdit-RaycastSomas
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${CMAKE_CURRENT_BINARY_DIR}/test_raycast_somas.xml
    ${TraceEdit_TEST_DIR}/baseline_raycast_somas.png
  )
  set_tests_properties(TraceEdit-RaycastSomas PROPERTIES
    DEPENDS TraceEdit-ExportCells
  )
  
  add_test(
    TraceEdit-DeleteTraces
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/cnic_macaque_pyramidal.swc
    ${TraceEdit_TEST_DIR}/test_delete_traces.xml
    ${TraceEdit_TEST_DIR}/baseline_delete_traces.png
  )
  set_tests_properties(TraceEdit-DeleteTraces PROPERTIES
    DEPENDS TraceEdit-RaycastSomas
  )
  
  add_test(
    TraceEdit-AddNewBranches
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_add_new_branches.xml
    ${TraceEdit_TEST_DIR}/baseline_add_new_branches.png
  )
  set_tests_properties(TraceEdit-AddNewBranches PROPERTIES
    DEPENDS TraceEdit-DeleteTraces
    PASS_REGULAR_EXPRESSION "Branching complete"
  )
  
  add_test(
    TraceEdit-CreateTraceBit
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_create_trace_bit.xml
    ${TraceEdit_TEST_DIR}/baseline_create_trace_bit.png
  )
  set_tests_properties(TraceEdit-CreateTraceBit PROPERTIES
    DEPENDS TraceEdit-AddNewBranches
  )
  
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_show_seed_points.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_show_seed_points.xml
    @ONLY)
  add_test(
    TraceEdit-ShowSeedPoints
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${CMAKE_CURRENT_BINARY_DIR}/test_show_seed_points.xml
    ${TraceEdit_TEST_DIR}/baseline_show_seed_points.png
  )
  set_tests_properties(TraceEdit-ShowSeedPoints PROPERTIES
    DEPENDS TraceEdit-AddNewBranches
  )
    
endif()
