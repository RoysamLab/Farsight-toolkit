ADD_EXECUTABLE(
  visualize
  ../Common/argparser.h
  ../Common/array.h
  ../Common/bag.h
  ../Common/boundingbox.cpp
  ../Common/boundingbox.h
  ../Common/camera.cpp
  ../Common/camera.h
  ../Common/draw_sphere.h
  ../Common/draw_sphere.cpp
  ../Common/edge.cpp
  ../Common/edge.h
  ../Common/face.cpp
  ../Common/face.h
  ../Common/gl_picking.cpp
  ../Common/glCanvas.cpp
  ../Common/glCanvas.h
  
  main.cpp
  ../Common/matrix.cpp
  ../Common/matrix.h
  ../Common/mesh.cpp
  ../Common/mesh.h
  ../Common/radiosity.cpp
  ../Common/radiosity.h
  ../Common/utils.h
  ../Common/vectors.h
  ../Common/vertex.h
  ../Common/vertex_parent.h
)
TARGET_LINK_LIBRARIES(visualize ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} vtksys)
INCLUDE_DIRECTORIES( ${GLUT_INCLUDE_DIR})



INSTALL( TARGETS visualize DESTINATION ${INSTALL_BIN_DIR} )

#Notes:set the Optimization option of the compiler(visual studio) as "Disabled/Od" other than "Maximize Speed/O2" before building