PROJECT(MDL)
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
SUBDIRS(ConnectTraces)
ADD_EXECUTABLE(AnisoDiffuse AnisoDiffuse.cpp)
ADD_EXECUTABLE(ConnCompntwFldfill ConnCompntwFldfill.cpp)
ADD_EXECUTABLE(GradientVecField GradientVecField.cpp)
ADD_EXECUTABLE(MinSpanTree distTransform.cpp MinSpanTree.cpp morphGraphPrune.cpp) 
ADD_EXECUTABLE(Integratedskel Integratedskel.cpp)
ADD_EXECUTABLE(volumeProcess volumeProcess.cpp)
ADD_EXECUTABLE(RefiningSkeleton RefiningSkeleton.cpp)
ADD_EXECUTABLE(BackboneExtract BackboneExtract.cpp)
TARGET_LINK_LIBRARIES(volumeProcess ITKIO ITKStatistics)

ADD_EXECUTABLE(AddSpinesToSmoothedBackbone AddSpinesToSmoothedBackbone.cxx)
TARGET_LINK_LIBRARIES(AddSpinesToSmoothedBackbone vtkCommon vtkIO vtkGraphics)

#INCLUDE_DIRECTORIES(${BSplineFitting_SOURCE_DIR})
ADD_EXECUTABLE(BSplineFitting BSplineFitting.cpp)
TARGET_LINK_LIBRARIES(BSplineFitting
ITKBasicFilters ITKCommon ITKIO ITKStatistics)

ADD_CUSTOM_COMMAND(
  TARGET MinSpanTree
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -D CMAKE_CFG_INTDIR:STRING=${CMAKE_CFG_INTDIR}
  ARGS -D CMAKE_EXECUTABLE_SUFFIX:STRING=${CMAKE_EXECUTABLE_SUFFIX}
  ARGS -D PROJECT_SOURCE_DIR:STRING=${PROJECT_SOURCE_DIR}
       -D PROJECT_BINARY_DIR:STRING=${PROJECT_BINARY_DIR}
       -P ${PROJECT_SOURCE_DIR}/ConfigureMDL.cmake)

ADD_CUSTOM_COMMAND(TARGET MinSpanTree POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  ${PROJECT_SOURCE_DIR}/ParseParameters.py
  ${PROJECT_BINARY_DIR}/ParseParameters.py)
ADD_CUSTOM_COMMAND(TARGET MinSpanTree POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  ${PROJECT_SOURCE_DIR}/MDLParameterDescriptions.xml
  ${PROJECT_BINARY_DIR}/MDLParameterDescriptions.xml)

INSTALL( TARGETS AnisoDiffuse DESTINATION bin/MDL )
INSTALL( TARGETS ConnCompntwFldfill DESTINATION bin/MDL )
INSTALL( TARGETS GradientVecField DESTINATION bin/MDL )
INSTALL( TARGETS MinSpanTree DESTINATION bin/MDL )
INSTALL( TARGETS Integratedskel DESTINATION bin/MDL )
INSTALL( TARGETS volumeProcess DESTINATION bin/MDL )
INSTALL( TARGETS BSplineFitting DESTINATION bin/MDL )
INSTALL( TARGETS RefiningSkeleton DESTINATION bin/MDL )
INSTALL( TARGETS BackboneExtract DESTINATION bin/MDL )
INSTALL( TARGETS AddSpinesToSmoothedBackbone DESTINATION bin/MDL )

