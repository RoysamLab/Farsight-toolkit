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
    ${TraceEdit_TEST_DIR}/test_reload.xml
    ${Farsight_BINARY_DIR}/Testing/Temporary/baseline_loadproject.png
    ${Farsight_BINARY_DIR}/Testing/Temporary/project_test.xml
  )
  set_tests_properties(TraceEdit-LoadProject PROPERTIES
    DEPENDS TraceEdit-SaveProject
  )
    
endif()
