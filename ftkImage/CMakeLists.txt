SET( FTKIMAGE_SRCS    
	vtkLSMReader.cxx
	ftkImage.cpp
	
)

SET( FTKIMAGE_HDRS
	vtkLSMReader.h
	vtkBXDProcessingWin32Header.h
	ftkImage.h
	ftkImage.txx
)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj" )

ADD_LIBRARY( ftkImage
             ${FTKIMAGE_SRCS} 
             ${FTKIMAGE_HDRS} 
)

TARGET_LINK_LIBRARIES( ftkImage ITKCommon ITKBasicFilters ITKIO vtksys vtkCommon vtkFiltering vtkIO vtkImaging )

#INSTALL( TARGETS ftkImage DESTINATION . )
                           
