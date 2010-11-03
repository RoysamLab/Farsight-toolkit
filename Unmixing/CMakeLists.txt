PROJECT(SpectralUnmixing)

cmake_minimum_required(VERSION 2.6)

SET(SRCS unmix_clustering.cpp)

ADD_EXECUTABLE(main ${SRCS})

TARGET_LINK_LIBRARIES(main ITKBasicFilters ITKCommon ITKIO itkvnl_inst ITKAlgorithms ITKgdcm ITKFEM ITKNumerics ITKStatistics ITKsys )

