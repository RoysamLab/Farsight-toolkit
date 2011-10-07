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
IF(WIN32)
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj" )
ENDIF(WIN32)

ADD_LIBRARY( ftkImage
             ${FTKIMAGE_SRCS} 
             ${FTKIMAGE_HDRS} 
)

TARGET_LINK_LIBRARIES( ftkImage vtksys vtkCommon vtkFiltering vtkIO vtkImaging ${ITK_LIBRARIES} ${OPENSLIDE_LIBRARY} )

#INSTALL( TARGETS ftkImage DESTINATION . )
                           
