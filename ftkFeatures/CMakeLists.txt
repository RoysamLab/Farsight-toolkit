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
	
)


ADD_LIBRARY( ftkFeatures
            ${FTKFEATURES_SRCS}
            ${FTKFEATURES_HDRS}
)



SUBDIRS(Qhull)


TARGET_LINK_LIBRARIES( ftkFeatures TinyXML ftkImage
                       vtksys vtkCommon vtkFiltering vtkIO
                       ${ITK_LIBRARIES} QHULL)
                       
IF (BUILD_ZERNIKE)
    TARGET_LINK_LIBRARIES( ftkFeatures zernike)
ENDIF (BUILD_ZERNIKE)
                      

