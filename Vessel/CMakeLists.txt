PROJECT(surface)

ADD_LIBRARY(surface farsight_surface_segmentation.cpp farsight_surface_segmentation.h parse_xml.cpp parse_xml.h)
TARGET_LINK_LIBRARIES(surface fscommon ITKCommon ITKIO ${LIBXML2_LIBRARIES} )
