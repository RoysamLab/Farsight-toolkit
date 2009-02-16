SET( cb_SRCS
  label2border.h
  itkLabelBorderImageFilter.h
  label2border.cpp
  itkLabelBorderImageFilter.txx
)

ADD_LIBRARY( cellLabelToBorders ${cb_SRCS} )
TARGET_LINK_LIBRARIES( cellLabelToBorders ITKBasicFilters ITKCommon ITKIO )
#INSTALL( TARGETS cellLabelToBorders DESTINATION . )
