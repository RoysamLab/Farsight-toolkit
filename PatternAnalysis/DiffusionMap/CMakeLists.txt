SET( DIFFUSION_SRCS
  DiffusionMap.cpp
)

SET( DIFFUSION_HDRS
  DiffusionMap.h
)


ADD_LIBRARY(diffusion_map ${DIFFUSION_SRCS} ${DIFFUSION_HDRS})

TARGET_LINK_LIBRARIES(diffusion_map ${VTK_LIBRARIES} ftkGraphs vnl vnl_algo vnl_io mbl)



