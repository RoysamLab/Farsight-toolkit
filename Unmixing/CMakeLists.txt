PROJECT(SpectralUnmixing)

cmake_minimum_required(VERSION 2.6)

SET(SRCS unmix_clustering.cpp)

ADD_EXECUTABLE(main ${SRCS})

TARGET_LINK_LIBRARIES(main ${ITK_LIBRARIES} )

