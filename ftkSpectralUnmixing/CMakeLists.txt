SET( SPECTRALUNMIXING_SRCS
  ftkSpectralUnmixing.cpp
)

SET( SPECTRALUNMIXING_HDRS
  ftkSpectralUnmixing.h
)


ADD_LIBRARY( ftkSpectralUnmixing
            ${SPECTRALUNMIXING_SRCS}
            ${SPECTRALUNMIXING_HDRS} 
)

TARGET_LINK_LIBRARIES( ftkSpectralUnmixing ${ITK_LIBRARIES} ftkImage)


