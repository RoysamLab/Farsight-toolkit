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


ADD_LIBRARY( ftkImage
             ${FTKIMAGE_SRCS} 
             ${FTKIMAGE_HDRS} 
)

TARGET_LINK_LIBRARIES( ftkImage vtksys vtkCommon vtkFiltering ITKCommon ITKBasicFilters ITKIO vtkIO vtkImaging )

#INSTALL( TARGETS ftkImage DESTINATION . )
                           
