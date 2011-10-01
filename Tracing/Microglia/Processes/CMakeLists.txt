add_executable( MicrogliaProcessTracer
  MPT_MAIN.cpp MicrogliaProcessTracer.cpp MicrogliaProcessTracer.h)

target_link_libraries(MicrogliaProcessTracer ${ITK_LIBRARIES})
