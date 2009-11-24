SET( FTKCOMMON_SRCS
  #ftkSegmentationResult.cpp
  ftkUtils.cpp
  ftkLabelImageToFeatures.cpp
  ftkIntrinsicFeatures.cpp
  ftkObject.cpp 
  ftkObjectAssociation.cpp
)

SET( FTKCOMMON_HDRS
  #ftkSegmentationResult.h
  ftkUtils.h
  ftkLabelImageToFeatures.h
  ftkLabelImageToFeatures.txx
  ftkIntrinsicFeatures.h
  ftkObject.h
  ftkObjectAssociation.h
  itkLabelGeometryImageFilter.h
  itkLabelGeometryImageFilter.txx 
)


ADD_LIBRARY( ftkCommon
            ${FTKCOMMON_SRCS} 
            ${FTKCOMMON_HDRS} 
)

TARGET_LINK_LIBRARIES( ftkCommon TinyXML
                       vtksys vtkCommon vtkFiltering
                       ITKCommon ITKBasicFilters ITKIO ITKStatistics vtkIO)

#INSTALL( TARGETS ftkCommon DESTINATION . )

