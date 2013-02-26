SET( MCLR_SRCS
  mclr.cpp
  mclr_SM.cpp
)

SET( MCLR_HDRS
  mclr.h
  mclr_SM.h
)

ADD_LIBRARY(MCLR ${MCLR_SRCS} ${MCLR_HDRS})
TARGET_LINK_LIBRARIES(MCLR
  ${VTK_LIBRARIES}
  vnl vnl_algo vnl_io mbl)
