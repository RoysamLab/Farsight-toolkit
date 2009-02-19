ADD_EXECUTABLE(centerline argparser.h array.h bag.h boundingbox.cpp boundingbox.h camera.cpp camera.h edge.cpp edge.h face.cpp face.h gl_picking.cpp glCanvas.cpp glCanvas.h main.cpp matrix.cpp matrix.h mesh.cpp mesh.h radiosity.cpp radisity.h utils.h vectors.h vertex.h vertex_parent.h)
ADD_EXECUTABLE(votestotiff votestotiff.cpp)

INSTALL(TARGETS centerline DESTINATION bin)
INSTALL(TARGETS votestotiff DESTINATION bin)

