ADD_EXECUTABLE(visualize argparser.h array.h bag.h boundingbox.cpp boundingbox.h camera.cpp camera.h edge.cpp edge.h face.cpp face.h gl_picking.cpp glCanvas.cpp glCanvas.h main.cpp matrix.cpp matrix.h mesh.cpp mesh.h radiosity.cpp radiosity.h utils.h vectors.h vertex.h vertex_parent.h)
INSTALL( TARGETS visualize DESTINATION bin)
