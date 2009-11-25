PROJECT(MDL)
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} ${Farsight_BINARY_DIR}/Tracing/MDL )
ADD_EXECUTABLE(AnisoDiffuse AnisoDiffuse.cpp)
ADD_EXECUTABLE(ConnCompntwFldfill ConnCompntwFldfill.cpp)
ADD_EXECUTABLE(GradientVecField GradientVecField.cpp)
ADD_EXECUTABLE(MinSpanTree distTransform.cpp MinSpanTree.cpp morphGraphPrune.cpp) 
ADD_EXECUTABLE(Integratedskel Integratedskel.cpp)
ADD_EXECUTABLE(volumeProcess volumeProcess.cpp)
TARGET_LINK_LIBRARIES(volumeProcess ITKIO ITKStatistics)
ADD_EXECUTABLE(RefiningSkeleton RefiningSkeleton.cpp)
ADD_EXECUTABLE(MDABasedSpineExtraction distTransform.cpp MDABasedSpineExtraction.cpp morphGraphPrune.cpp)
ADD_EXECUTABLE(BackboneExtract BackboneExtract.cpp)

#ADD_EXECUTABLE(AddSpinesToSmoothedBackbone AddSpinesToSmoothedBackbone.cxx)
#TARGET_LINK_LIBRARIES(AddSpinesToSmoothedBackbone vtkCommon vtkIO vtkGraphics)

ADD_EXECUTABLE(BSplineFitting BSplineFitting.cpp)
TARGET_LINK_LIBRARIES(BSplineFitting
  ITKBasicFilters ITKCommon ITKIO ITKStatistics)

QT4_WRAP_UI( MY_UI_HDRS MDLGUI.ui )
QT4_WRAP_CPP( GUI_MOC_SRCS MDLGUI.h )
SET(MDLGUI_SRCS
  MDLGUI.h MDLGUI.cxx 
  MDLGUI-main.cxx ${GUI_MOC_SRCS} ${MY_UI_HDRS})
ADD_EXECUTABLE(MDLGUI ${MDLGUI_SRCS} )
TARGET_LINK_LIBRARIES(MDLGUI ${QT_LIBRARIES})

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

INSTALL( TARGETS AnisoDiffuse DESTINATION bin/Tracing/MDL )
INSTALL( TARGETS ConnCompntwFldfill DESTINATION bin/Tracing/MDL )
INSTALL( TARGETS GradientVecField DESTINATION bin/Tracing/MDL )
INSTALL( TARGETS MinSpanTree DESTINATION bin/Tracing/MDL )
INSTALL( TARGETS Integratedskel DESTINATION bin/Tracing/MDL )
INSTALL( TARGETS volumeProcess DESTINATION bin/Tracing/MDL )
INSTALL( TARGETS BSplineFitting DESTINATION bin/Tracing/MDL )
INSTALL( TARGETS RefiningSkeleton DESTINATION bin/Tracing/MDL )
INSTALL( TARGETS MDLGUI DESTINATION bin/Tracing/MDL )
#INSTALL( TARGETS AddSpinesToSmoothedBackbone DESTINATION bin/Tracing/MDL )
INSTALL( TARGETS MDABasedSpineExtraction DESTINATION bin/Tracing/MDL)
