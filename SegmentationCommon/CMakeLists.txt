SET( SEGMENTATIONCOMMON_SRCS
  #ftkSegmentationResult.cpp
  ftkIntrinsicFeatures.cpp
  ftkObject.cpp 
  ftkObjectAssociation.cpp
)

SET( SEGMENTATIONCOMMON_HDRS
  #ftkSegmentationResult.h
  ftkLabelImageToFeatures.h
  ftkLabelImageToFeatures.txx
  ftkIntrinsicFeatures.h
  ftkObject.h
  ftkObjectAssociation.h
  itkLabelGeometryImageFilter.h
  itkLabelGeometryImageFilter.txx 
)


ADD_LIBRARY( SegmentationCommon
            ${SEGMENTATIONCOMMON_SRCS} 
            ${SEGMENTATIONCOMMON_HDRS} 
)

TARGET_LINK_LIBRARIES( SegmentationCommon TinyXML
                       vtksys vtkCommon vtkFiltering
                       ITKCommon ITKBasicFilters ITKIO ITKStatistics vtkIO)

#INSTALL( TARGETS SegmentationCommon DESTINATION . )

