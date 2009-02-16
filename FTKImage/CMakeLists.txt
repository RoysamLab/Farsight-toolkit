SET( FTKIMAGE_SRCS    
	vtkLSMReader.cxx
	vtkKWImage.cxx
	vtkKWImageIO.cxx
	ftkImage.cpp	
)

SET( FTKIMAGE_HDRS
	vtkLSMReader.h
	vtkBXDProcessingWin32Header.h
	vtkKWImage.h
	vtkKWImageIO.h
	ftkImage.h
	ftkImage.txx
)


ADD_LIBRARY( FTKImage 
             ${FTKIMAGE_SRCS} 
             ${FTKIMAGE_HDRS} 
)

TARGET_LINK_LIBRARIES( FTKImage vtksys vtkCommon vtkFiltering 
				ITKCommon ITKBasicFilters ITKIO vtkIO vtkImaging )

#INSTALL( TARGETS FTKImage DESTINATION . )
                           
