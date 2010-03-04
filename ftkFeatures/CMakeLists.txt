SET( FTKFEATURES_SRCS
  ftkLabelImageToFeatures.cpp
  ftkIntrinsicFeatures.cpp
  ftkObject.cpp 
  ftkObjectAssociation.cpp
)

SET( FTKFEATURES_HDRS
  ftkLabelImageToFeatures.h
  ftkLabelImageToFeatures.txx
  ftkIntrinsicFeatures.h
  ftkObject.h
  ftkObjectAssociation.h
  itkLabelGeometryImageFilter.h
  itkLabelGeometryImageFilter.txx 
)

ADD_LIBRARY( ftkFeatures
            ${FTKFEATURES_SRCS}
            ${FTKFEATURES_HDRS}
)

TARGET_LINK_LIBRARIES( ftkFeatures TinyXML ftkImage
                       vtksys vtkCommon vtkFiltering vtkIO
                       ITKCommon ITKBasicFilters ITKIO ITKStatistics )

