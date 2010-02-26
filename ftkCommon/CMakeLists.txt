
ADD_LIBRARY( TraceProject TraceProjectManager.h TraceProjectManager.cxx )
TARGET_LINK_LIBRARIES( TraceProject TinyXML )


SET( FTKCOMMON_SRCS
  #ftkSegmentationResult.cpp
  ftkUtils.cpp
  ftkLabelImageToFeatures.cpp
  ftkIntrinsicFeatures.cpp
  ftkObject.cpp 
  ftkObjectAssociation.cpp
  ftkPreprocess.cpp
  ftkMultipleImageHandler.cpp
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
  ftkPreprocess.h
  ftkMultipleImageHandler.h
)

ADD_LIBRARY( ftkCommon
            ${FTKCOMMON_SRCS} 
            ${FTKCOMMON_HDRS} 
)

TARGET_LINK_LIBRARIES( ftkCommon TinyXML ftkImage TraceProject
                       vtksys vtkCommon vtkFiltering vtkIO
                       ITKCommon ITKBasicFilters ITKIO ITKStatistics )

ADD_EXECUTABLE( block main.cpp )
TARGET_LINK_LIBRARIES( block ftkCommon )

#INSTALL( TARGETS ftkCommon DESTINATION . )

