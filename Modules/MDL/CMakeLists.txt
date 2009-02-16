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
ADD_CUSTOM_COMMAND(
  TARGET MinSpanTree
  PRE_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PROJECT_SOURCE_DIR}/Parameters.xml
               ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/Parameters.xml)
ADD_CUSTOM_COMMAND(
  TARGET MinSpanTree
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -D CMAKE_CFG_INTDIR:STRING=${CMAKE_CFG_INTDIR}
  ARGS -D CMAKE_EXECUTABLE_SUFFIX:STRING=${CMAKE_EXECUTABLE_SUFFIX}
  ARGS -D PROJECT_SOURCE_DIR:STRING=${PROJECT_SOURCE_DIR}
       -D PROJECT_BINARY_DIR:STRING=${PROJECT_BINARY_DIR}
       -P ${PROJECT_SOURCE_DIR}/ConfigureMDL.cmake)
