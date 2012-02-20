cmake_minimum_required(VERSION 2.8)
PROJECT(TraceEdit)

QT4_WRAP_CPP( GUI_MOC_SRCS TraceView3D.h MergeModel.h TraceModel.h CellTraceModel.h cellexport.h screenshot.h)
  
SET( TRACE_SRCS
  TraceBit.h        TraceBit.cxx
  TraceLine.h       TraceLine.cxx
  TraceGap.h        TraceGap.cxx
  TraceObject.h     TraceObject.cxx
  TraceModel.h      TraceModel.cxx
  MergeModel.h      MergeModel.cxx
  branchPT.h        BranchPT.cxx
  ImageActors.h     ImageActors.cxx
  CellTrace.h       CellTrace.cxx
  CellTraceModel.h  CellTraceModel.cxx
  vtkPlotEdges.h    vtkPlotEdges.cxx 
  TraceView3D.h     TraceView3D.cxx 
  cellexport.h	    cellexport.cxx
  screenshot.h	    screenshot.cxx
  GridActors.h		GridActors.cxx
  VolumeOfInterest.h VolumeOfInterest.cxx
  ${GUI_MOC_SRCS} ) 

SET( TRACE_LIBS
    tiff
    ${ITK_LIBRARIES}
    vtkRendering
    vtkIO
    vtkHybrid
    vtkWidgets
    vtkCommon
    vtksys
    vtkFiltering
    vtkVolumeRendering
    QVTK
    ftkGUI
    ftkCommon
    ftkGraphs
    MCLR
    ${QT_LIBRARIES} )

IF (BUILD_SPD)	
  LIST( APPEND TRACE_LIBS SPD ) 
ENDIF (BUILD_SPD)

IF(BUILD_CLUSCLUS)
  LIST( APPEND TRACE_LIBS ClusClus_lib )
ENDIF (BUILD_CLUSCLUS)

IF(QtTestingFound)
  LINK_DIRECTORIES( ${QtTesting_LIBRARY_DIR} )
  ADD_DEFINITIONS(-DUSE_QT_TESTING)
ENDIF()

OPTION( USE_GPUREN "Use GPU for Rendering" OFF )
IF(USE_GPUREN)
	ADD_DEFINITIONS(-DUSE_GPUREN)
ENDIF(USE_GPUREN)

ADD_LIBRARY( Trace STATIC ${TRACE_SRCS} )
TARGET_LINK_LIBRARIES( Trace ${TRACE_LIBS} )

ADD_EXECUTABLE(${PROJECT_NAME} ${EXE_TYPE} vtkTraceViewerDriver.cxx )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} Trace)

if(BUILD_TESTING)
  add_subdirectory(Testing)
endif(BUILD_TESTING)

if(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
   set(plugin_path ${FARSIGHT_BUNDLE_LOCATION}/${PROJECT_NAME}.app/Contents/plugins/imageformats)
   set(resources_path ${FARSIGHT_BUNDLE_LOCATION}/${PROJECT_NAME}.app/Contents/Resources)
else(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION bin )
   set(plugin_path "${CMAKE_INSTALL_PREFIX}/bin")
   set(resources_path "${CMAKE_INSTALL_PREFIX}/bin")
endif(APPLE)

#install qt.conf & Qt plugins loaded by TraceEdit at runtime
install(FILES ${QT_QGIF_PLUGIN_RELEASE} DESTINATION ${plugin_path})
install(FILES ${QT_QICO_PLUGIN_RELEASE} DESTINATION ${plugin_path})
install(FILES ${QT_QJPEG_PLUGIN_RELEASE} DESTINATION ${plugin_path})
install(FILES ${QT_QMNG_PLUGIN_RELEASE} DESTINATION ${plugin_path})
install(FILES ${QT_QSVG_PLUGIN_RELEASE} DESTINATION ${plugin_path})
install(FILES ${QT_QTIFF_PLUGIN_RELEASE} DESTINATION ${plugin_path})
install(FILES ${Farsight_SOURCE_DIR}/InstallRules/qt.conf DESTINATION ${resources_path})

