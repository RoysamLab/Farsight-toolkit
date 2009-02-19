
ADD_EXECUTABLE(centerline
  argparser.h
  array.h
  bag.h
  boundingbox.cpp
  boundingbox.h
  camera.cpp
  camera.h
  draw_sphere.h
  draw_sphere.cpp
  edge.cpp
  edge.h
  face.cpp  
  face.h
  gl_picking.cpp
  glCanvas.cpp
  glCanvas.h
  main.cpp
  matrix.cpp
  matrix.h
  mesh.cpp
  mesh.h
  radiosity.cpp
  radiosity.h
  utils.h
  vectors.h
  vertex.h
  vertex_parent.h)

TARGET_LINK_LIBRARIES(centerline ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

ADD_EXECUTABLE(votestotiff votestotiff.cpp)
TARGET_LINK_LIBRARIES(votestotiff ITKCommon ITKIO)

INSTALL(TARGETS centerline DESTINATION bin)
INSTALL(TARGETS votestotiff DESTINATION bin)

