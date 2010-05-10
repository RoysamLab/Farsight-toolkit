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

ADD_LIBRARY( ftkPreprocess
	ftkPreprocess2.h
	ftkPreprocess2.cpp
	itkBinaryThinningImageFilter3D.h
)
TARGET_LINK_LIBRARIES( ftkPreprocess ITKCommon ITKBasicFilters ITKIO ITKStatistics TinyXML)

ADD_EXECUTABLE( block main.cpp )
TARGET_LINK_LIBRARIES( block ftkCommon )

ADD_EXECUTABLE( prep preprocess_main.cpp )
TARGET_LINK_LIBRARIES( prep ftkPreprocess )

#INSTALL( TARGETS ftkCommon DESTINATION . )

