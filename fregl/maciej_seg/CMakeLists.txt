
SET( MS_sources
  pix_t.cpp 
  wtd_lg.cpp 
  filter.h                filter.cpp 
  cell_mgr.h              cell_mgr.cpp
  cell.h 
  pix_t.h 
  wtd_lg.h 
  compare.h
  maciejSegmentation.h    maciejSegmentation.cpp
  xml_util.h              xml_util.cpp
  rich_cell.h             rich_cell.cpp
  farsight_maciejseg.h    farsight_maciejseg.cpp  
  #nucleus_feature_util.h  nucleus_feature_util.cpp 
  fnsl3d_filter.h         fnsl3d_filter.cxx
  nuclear_IDL_util.h      nuclear_IDL_util.cpp
)

AUX_SOURCE_DIRECTORY( Templates MS_sources)

ADD_LIBRARY( maciejseg ${MS_sources} )

TARGET_LINK_LIBRARIES( maciejseg fscommon rsdl vbl  
${LIBXML2_LIBRARIES}  ${ITK_LIBRARIES} )
