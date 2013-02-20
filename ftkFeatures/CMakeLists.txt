if(BUILD_ZERNIKE)
  if(MSVC)
    add_definitions( -DZERNIKE )
  endif(MSVC)
  include_directories( ${Farsight_SOURCE_DIR}/Zernike )
endif(BUILD_ZERNIKE)

set( FTKFEATURES_SRCS
  ftkLabelImageToFeatures.cpp
  ftkIntrinsicFeatures.cpp
  ftkObject.cpp
  ftkObjectAssociation.cpp
  )

set( FTKFEATURES_HDRS
  ftkLabelImageToFeatures.h
  ftkLabelImageToFeatures.txx
  ftkIntrinsicFeatures.h
  ftkObject.h
  ftkObjectAssociation.h
  )


add_library( ftkFeatures
  ${FTKFEATURES_SRCS}
  ${FTKFEATURES_HDRS}
  )
add_subdirectory(Qhull)
target_link_libraries( ftkFeatures
  TinyXML ftkImage
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  QHULL
  )
if( BUILD_ZERNIKE )
  target_link_libraries( ftkFeatures zernike )
endif()
