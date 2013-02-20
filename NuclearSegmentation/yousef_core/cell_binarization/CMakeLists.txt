SET( cb_SRCS
  cell_binarization.h
  cell_binarization.cpp    
  adaptive_binarization.h
  adaptive_binarization.cpp
  new_graph.h
  block_b.h   
)

ADD_LIBRARY( cell_binarization ${cb_SRCS} )
TARGET_LINK_LIBRARIES ( cell_binarization        ${ITK_LIBRARIES} )

#INSTALL( TARGETS cell_binarization DESTINATION ${INSTALL_BIN_DIR} )
