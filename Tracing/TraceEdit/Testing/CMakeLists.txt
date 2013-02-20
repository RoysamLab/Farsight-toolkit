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
    if(UNIX)
      set(AppName "TraceEdit-real")
    else()
      set(AppName "TraceEdit")
    endif()
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
    DEPENDS TraceEdit-CreateTraceBit
    PASS_REGULAR_EXPRESSION "seed points rendered"
  )
  
  add_test(
    TraceEdit-CellAnalysis
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_cell_analysis.xml
    ${TraceEdit_TEST_DIR}/baseline_add_new_branches.png
  )
  set_tests_properties(TraceEdit-CellAnalysis PROPERTIES
    DEPENDS TraceEdit-ShowSeedPoints
  )
  
  add_test(
    TraceEdit-RenderModes
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    ${TraceEdit_TEST_DIR}/test_render_modes.xml
    ${TraceEdit_TEST_DIR}/baseline_render_modes.png
  )
  set_tests_properties(TraceEdit-RenderModes PROPERTIES
    DEPENDS TraceEdit-CellAnalysis
  )

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_save_to_file.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_save_to_file.xml
    @ONLY)
  add_test(
    TraceEdit-SaveToFile
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${CMAKE_CURRENT_BINARY_DIR}/test_save_to_file.xml
    ${TraceEdit_TEST_DIR}/baseline_add_new_branches.png
  )
  set_tests_properties(TraceEdit-SaveToFile PROPERTIES
    DEPENDS TraceEdit-RenderModes
    PASS_REGULAR_EXPRESSION "File saved as"
  )
  
  add_test(
    TraceEdit-ShowPoints
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_show_points.xml
    ${TraceEdit_TEST_DIR}/baseline_show_points.png
  )
  set_tests_properties(TraceEdit-ShowPoints PROPERTIES
    DEPENDS TraceEdit-SaveToFile
  )
  
  add_test(
    TraceEdit-SetRoots
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_set_roots.xml
    ${TraceEdit_TEST_DIR}/baseline_add_new_branches.png
  )
  set_tests_properties(TraceEdit-SetRoots PROPERTIES
    DEPENDS TraceEdit-RenderModes
    PASS_REGULAR_EXPRESSION "set roots"
  )

  add_test(
    TraceEdit-ListSelections
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_list_selections.xml
    ${TraceEdit_TEST_DIR}/baseline_add_new_branches.png
  )
  set_tests_properties(TraceEdit-ListSelections PROPERTIES
    DEPENDS TraceEdit-SetRoots
    PASS_REGULAR_EXPRESSION "65 lines are selected"
  )
  
  add_test(
    TraceEdit-RemoveImageActors
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    ${TraceEdit_TEST_DIR}/test_remove_image_actors.xml
    ${TraceEdit_TEST_DIR}/baseline_delete_traces.png
  )
  set_tests_properties(TraceEdit-RemoveImageActors PROPERTIES
    DEPENDS TraceEdit-SetRoots
  )
  
  add_test(
    TraceEdit-BreakBranch
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_break_branch.xml
    ${TraceEdit_TEST_DIR}/baseline_break_branch.png
  )
  set_tests_properties(TraceEdit-BreakBranch PROPERTIES
    DEPENDS TraceEdit-RemoveImageActors
  )
  
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_voi.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_voi.xml
    @ONLY)
  add_test(
    TraceEdit-VOI
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${CMAKE_CURRENT_BINARY_DIR}/test_voi.xml
    ${TraceEdit_TEST_DIR}/baseline_voi.png
  )
  set_tests_properties(TraceEdit-VOI PROPERTIES
    DEPENDS TraceEdit-BreakBranch
  )
  
  add_test(
    TraceEdit-SetTracetype
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_set_trace_type.xml
    ${TraceEdit_TEST_DIR}/baseline_set_trace_type.png
  )
  set_tests_properties(TraceEdit-SetTracetype PROPERTIES
    DEPENDS TraceEdit-VOI
  )
  
  add_test(
    TraceEdit-ShowStatistics
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_show_statistics.xml
    ${TraceEdit_TEST_DIR}/baseline_add_new_branches.png
  )
  set_tests_properties(TraceEdit-ShowStatistics PROPERTIES
    DEPENDS TraceEdit-VOI
  )
  
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_vtk_traces.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_vtk_traces.xml
    @ONLY)
  add_test(
    TraceEdit-VTKTraces
    ${TraceEdit_EXE}
    ${CMAKE_CURRENT_BINARY_DIR}/test_vtk_traces.xml
    ${FARSIGHT_DATA_ROOT}/TraceData/trach11a.vtk
    ${TraceEdit_TEST_DIR}/baseline_vtk_traces.png
  )
  set_tests_properties(TraceEdit-VTKTraces PROPERTIES
    DEPENDS TraceEdit-ShowStatistics
    PASS_REGULAR_EXPRESSION "trach11a.vtk"
  )
  
  add_test(
    TraceEdit-Gridlines-XY
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_grid_xy.xml
    ${TraceEdit_TEST_DIR}/baseline_grid_xy.png
  )
  set_tests_properties(TraceEdit-Gridlines-XY PROPERTIES
    DEPENDS TraceEdit-VTKTraces
  )
  
  add_test(
    TraceEdit-Gridlines-XZ
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_grid_xz.xml
    ${TraceEdit_TEST_DIR}/baseline_grid_xz.png
  )
  set_tests_properties(TraceEdit-Gridlines-XZ PROPERTIES
    DEPENDS TraceEdit-Gridlines-XY
  )

  add_test(
    TraceEdit-Gridlines-YZ
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop-trace.swc
    ${TraceEdit_TEST_DIR}/test_grid_yz.xml
    ${TraceEdit_TEST_DIR}/baseline_grid_yz.png
  )
  set_tests_properties(TraceEdit-Gridlines-YZ PROPERTIES
    DEPENDS TraceEdit-Gridlines-XZ
  )
  
  add_test(
    TraceEdit-ColorByTrees
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/9cells.swc
    ${TraceEdit_TEST_DIR}/test_color_by_trees.xml
    ${TraceEdit_TEST_DIR}/baseline_color_by_trees.png
  )
  set_tests_properties(TraceEdit-ColorByTrees PROPERTIES
    DEPENDS TraceEdit-Gridlines-YZ
  )
    
endif()
