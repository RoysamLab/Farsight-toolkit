if(QtTestingFound)
  if(APPLE)
    set(NucleusEditor_EXE ${Farsight_BINARY_DIR}/exe/NucleusEditor.app/Contents/MacOS/NucleusEditor)
    set(AppName "Farsight")
  else()
    set(NucleusEditor_EXE ${Farsight_BINARY_DIR}/exe/NucleusEditor)
    if(UNIX)
      set(AppName "NucleusEditor-real")
    else()
      set(AppName "NucleusEditor")
    endif()
  endif()

  set(NucleusEditor_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

  add_test(
    NucleusEditor-ROIStatistics
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${NucleusEditor_TEST_DIR}/ROI_statistics.xml
  )
  set_tests_properties(NucleusEditor-ROIStatistics PROPERTIES
    PASS_REGULAR_EXPRESSION "variance:"
  )
  
  add_test(
    NucleusEditor-SliceViewer
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${NucleusEditor_TEST_DIR}/slice_viewer.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/slice_viewer_baseline.png
  )
  set_tests_properties(NucleusEditor-SliceViewer PROPERTIES
    DEPENDS NucleusEditor-ROIStatistics
  )
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/segment_nuclei.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/segment_nuclei.xml
                 @ONLY)
  add_test(
    NucleusEditor-SegmentNuclei
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/segment_nuclei.xml
  )
  set_tests_properties(NucleusEditor-SegmentNuclei PROPERTIES
    PASS_REGULAR_EXPRESSION "Cell Binarized.. with 11 connected components"
    DEPENDS NucleusEditor-SliceViewer
  )

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/merge_segments.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/merge_segments.xml
                 @ONLY)
  add_test(
    NucleusEditor-MergeSegments
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/merge_segments.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/merge_segments_baseline.png
  )
  set_tests_properties(NucleusEditor-MergeSegments PROPERTIES
    DEPENDS NucleusEditor-SegmentNuclei
  )
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/delete_segments.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/delete_segments.xml
                 @ONLY)
  add_test(
    NucleusEditor-DeleteSegments
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/delete_segments.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/delete_segments_baseline.png
  )
  set_tests_properties(NucleusEditor-DeleteSegments PROPERTIES
    DEPENDS NucleusEditor-MergeSegments
  )
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/show_centroids.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/show_centroids.xml
                 @ONLY)
  add_test(
    NucleusEditor-ShowCentroids
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/show_centroids.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/show_centroids_baseline.png
  )
  set_tests_properties(NucleusEditor-ShowCentroids PROPERTIES
    DEPENDS NucleusEditor-DeleteSegments
  )

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/show_object_IDs.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/show_object_IDs.xml
                 @ONLY)
  add_test(
    NucleusEditor-ShowObjectIDs
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/show_object_IDs.xml
  )
  set_tests_properties(NucleusEditor-ShowObjectIDs PROPERTIES
    DEPENDS NucleusEditor-ShowCentroids
		PASS_REGULAR_EXPRESSION "toggleIDs"
  )

  add_test(
    NucleusEditor-About
    ${NucleusEditor_EXE}
    --test
    ${NucleusEditor_TEST_DIR}/about.xml
  )
  set_tests_properties(NucleusEditor-About PROPERTIES
    DEPENDS NucleusEditor-ShowObjectIDs
		PASS_REGULAR_EXPRESSION "About FARSIGHT"
  )
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/exit.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/exit.xml
                 @ONLY)
  add_test(
    NucleusEditor-Exit
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/exit.xml
  )
  set_tests_properties(NucleusEditor-Exit PROPERTIES
    DEPENDS NucleusEditor-About
		PASS_REGULAR_EXPRESSION "closing NucleusEditor"
  )

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/detect_outliers.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/detect_outliers.xml
                 @ONLY)
  add_test(
    NucleusEditor-DetectOutliers
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/detect_outliers.xml
  )
  set_tests_properties(NucleusEditor-DetectOutliers PROPERTIES
    DEPENDS NucleusEditor-Exit
		PASS_REGULAR_EXPRESSION "optimization finished, #iter = 0"
  )
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/save_project.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/save_project.xml
                 @ONLY)
  add_test(
    NucleusEditor-SaveProject
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/save_project.xml
  )
  set_tests_properties(NucleusEditor-SaveProject PROPERTIES
    DEPENDS NucleusEditor-DetectOutliers
		PASS_REGULAR_EXPRESSION "done saving"
  )

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/load_project.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/load_project.xml
                 @ONLY)
  add_test(
    NucleusEditor-LoadProject
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/load_project.xml
  )
  set_tests_properties(NucleusEditor-LoadProject PROPERTIES
    DEPENDS NucleusEditor-SaveProject
		PASS_REGULAR_EXPRESSION "project loaded"
  )

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/update_database.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/update_database.xml
                 @ONLY)
  add_test(
    NucleusEditor-UpdateDatabase
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/update_database.xml
  )
  set_tests_properties(NucleusEditor-UpdateDatabase PROPERTIES
    DEPENDS NucleusEditor-LoadProject
		PASS_REGULAR_EXPRESSION "The image ID on the database is: 1"
  )

  add_test(
    NucleusEditor-StartActiveLearning
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${NucleusEditor_TEST_DIR}/start_active_learning.xml
  )
  set_tests_properties(NucleusEditor-StartActiveLearning PROPERTIES
    DEPENDS NucleusEditor-UpdateDatabase
		PASS_REGULAR_EXPRESSION "Max Info-0.0357334"
    TIMEOUT 1800) #30 minutes
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/unmix_channels.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/unmix_channels.xml
                 @ONLY)
  add_test(
    NucleusEditor-UnmixChannels
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/unmix_channels.xml
  )
  set_tests_properties(NucleusEditor-UnmixChannels PROPERTIES
    DEPENDS NucleusEditor-StartActiveLearning
		PASS_REGULAR_EXPRESSION "done with unmixing"
    TIMEOUT 1800) #30 minutes
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/clustering_heatmap.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/clustering_heatmap.xml
                 @ONLY)
  add_test(
    NucleusEditor-ClusteringHeatmap
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/clustering_heatmap.xml
  )
  set_tests_properties(NucleusEditor-ClusteringHeatmap PROPERTIES
    DEPENDS NucleusEditor-UnmixChannels
		PASS_REGULAR_EXPRESSION "get selected"
    TIMEOUT 1800) #30 minutes
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fill_cells.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/fill_cells.xml
                 @ONLY)
  add_test(
    NucleusEditor-FillCells
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/fill_cells.xml
  )
  set_tests_properties(NucleusEditor-FillCells PROPERTIES
    DEPENDS NucleusEditor-ClusteringHeatmap
		PASS_REGULAR_EXPRESSION "cells filled")

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/split_cell.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/split_cell.xml
                 @ONLY)
  add_test(
    NucleusEditor-SplitCell
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/split_cell.xml
  )
  set_tests_properties(NucleusEditor-SplitCell PROPERTIES
    DEPENDS NucleusEditor-FillCells
		PASS_REGULAR_EXPRESSION "cell split")
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/exclude_margins.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/exclude_margins.xml
                 @ONLY)
  add_test(
    NucleusEditor-ExcludeMargins
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/exclude_margins.xml
  )
  set_tests_properties(NucleusEditor-ExcludeMargins PROPERTIES
    DEPENDS NucleusEditor-SplitCell
		PASS_REGULAR_EXPRESSION "margins excluded")
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/classify.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/classify.xml
                 @ONLY)
  add_test(
    NucleusEditor-Classify
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/classify.xml
  )
  set_tests_properties(NucleusEditor-Classify PROPERTIES
    DEPENDS NucleusEditor-ExcludeMargins
    PASS_REGULAR_EXPRESSION "KPLS classification complete")
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/change_class.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/change_class.xml
                 @ONLY)
  add_test(
    NucleusEditor-ChangeClass
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/change_class.xml
  )
  set_tests_properties(NucleusEditor-ChangeClass PROPERTIES
    DEPENDS NucleusEditor-Classify
    PASS_REGULAR_EXPRESSION "class changed for 3 objects.")
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/KNN.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/KNN.xml
                 @ONLY)
  add_test(
    NucleusEditor-KNN
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/KNN.xml
  )
  set_tests_properties(NucleusEditor-KNN PROPERTIES
    PASS_REGULAR_EXPRESSION "5 neighbors potentially connected"
    DEPENDS NucleusEditor-Classify)
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/query_radius.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/query_radius.xml
                 @ONLY)
  add_test(
    NucleusEditor-QueryRadius
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/query_radius.xml
  )
  set_tests_properties(NucleusEditor-QueryRadius PROPERTIES
    PASS_REGULAR_EXPRESSION "3 neighbors in radius"
    DEPENDS NucleusEditor-KNN)
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/split_along_z.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/split_along_z.xml
                 @ONLY)
  add_test(
    NucleusEditor-SplitAlongZ
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/split_along_z.xml
  )
  set_tests_properties(NucleusEditor-SplitAlongZ PROPERTIES
    PASS_REGULAR_EXPRESSION "split along Z"
    DEPENDS NucleusEditor-QueryRadius)
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/batch_split.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/batch_split.xml
                 @ONLY)
  add_test(
    NucleusEditor-BatchSplit
    ${NucleusEditor_EXE}
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/batch_split.xml
  )
  set_tests_properties(NucleusEditor-BatchSplit PROPERTIES
    PASS_REGULAR_EXPRESSION "batch split"
    DEPENDS NucleusEditor-SplitAlongZ)

endif()

