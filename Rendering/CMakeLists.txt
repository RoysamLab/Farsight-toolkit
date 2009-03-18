# set project's name
PROJECT( Rendering )

  FIND_PACKAGE(VTK)
  IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
  ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR
    "Cannot build Farsight without VTK.  Please set VTK_DIR.")
  ENDIF(VTK_FOUND)

  FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
    "Cannot build Farsight without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


ADD_EXECUTABLE(render main.cpp)
TARGET_LINK_LIBRARIES(render vtkRendering vtkCommon vtkIO vtksys vtkFiltering vtkHybrid)
TARGET_LINK_LIBRARIES(render ITKBasicFilters ITKCommon ITKIO)
INSTALL(TARGETS render DESTINATION bin)
