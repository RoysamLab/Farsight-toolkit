PROJECT(MDL)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ADD_EXECUTABLE(AnisoDiffuse AnisoDiffuse.cpp)
ADD_EXECUTABLE(ConnCompntwFldfill ConnCompntwFldfill.cpp)
ADD_EXECUTABLE(ConnectComponents ConnectComponents.cpp)
ADD_EXECUTABLE(floodfill floodfill.cpp)
ADD_EXECUTABLE(Gen_VectorField Gen_VectorField.cpp)
ADD_EXECUTABLE(GradientVecField GradientVecField.cpp)
ADD_EXECUTABLE(MinSpanTree distTransform.cpp MinSpanTree.cpp morphGraphPrune.cpp) 
ADD_EXECUTABLE(skel_Extrvalley3D skel_Extrvalley3D.cpp)
ADD_EXECUTABLE(skel_streamline skel_streamline.cpp)
ADD_EXECUTABLE(volumeProcess volumeProcess.cpp)
#ADD_CUSTOM_COMMAND(
#  TARGET MinSpanTree
#  PRE_BUILD
#  COMMAND ${CMAKE_COMMAND}
#  ARGS -E copy ${Farsight_SOURCE_DIR}/python/XML/MDLParameterDescriptions.xml
#               ${Farsight_BINARY_DIR}/python/XML/MDLParameterDescriptions.xml)
ADD_CUSTOM_COMMAND(
  TARGET MinSpanTree
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -D CMAKE_CFG_INTDIR:STRING=${CMAKE_CFG_INTDIR}
  ARGS -D CMAKE_EXECUTABLE_SUFFIX:STRING=${CMAKE_EXECUTABLE_SUFFIX}
  ARGS -D PROJECT_SOURCE_DIR:STRING=${PROJECT_SOURCE_DIR}
       -D Farsight_SOURCE_DIR:STRING=${Farsight_SOURCE_DIR}
       -P ${PROJECT_SOURCE_DIR}/ConfigureMDL.cmake)

INSTALL( TARGETS AnisoDiffuse DESTINATION bin/MDL )
INSTALL( TARGETS ConnCompntwFldfill DESTINATION bin/MDL )
INSTALL( TARGETS ConnectComponents DESTINATION bin/MDL )
INSTALL( TARGETS floodfill DESTINATION bin/MDL )
INSTALL( TARGETS Gen_VectorField DESTINATION bin/MDL )
INSTALL( TARGETS GradientVecField DESTINATION bin/MDL )
INSTALL( TARGETS MinSpanTree DESTINATION bin/MDL )
INSTALL( TARGETS skel_Extrvalley3D DESTINATION bin/MDL )
INSTALL( TARGETS skel_streamline DESTINATION bin/MDL )
INSTALL( TARGETS volumeProcess DESTINATION bin/MDL )
