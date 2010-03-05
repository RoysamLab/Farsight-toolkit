SET( FTKCOMMON_SRCS
  ftkUtils.cpp
  ftkPreprocess.cpp
  ftkMultipleImageHandler.cpp
  ftkParameters.cpp
  ftkProjectManager.cxx
)

SET( FTKCOMMON_HDRS
  ftkUtils.h
  ftkPreprocess.h
  ftkMultipleImageHandler.h
  ftkParameters.h
  ftkProjectManager.h
)

ADD_LIBRARY( ftkCommon
            ${FTKCOMMON_SRCS} 
            ${FTKCOMMON_HDRS} 
)

TARGET_LINK_LIBRARIES( ftkCommon TinyXML ftkImage
                       vtksys vtkCommon vtkFiltering vtkIO
                       ITKCommon ITKBasicFilters ITKIO ITKStatistics )

ADD_EXECUTABLE( block main.cpp )
TARGET_LINK_LIBRARIES( block ftkCommon )

#INSTALL( TARGETS ftkCommon DESTINATION . )

