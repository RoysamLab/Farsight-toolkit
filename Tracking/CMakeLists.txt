PROJECT(Tracker)

cmake_minimum_required(VERSION 2.0)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)



SET(SRCS helpers.h itkBinaryThinningImageFilter3D.h ftkTrackFeatures.h itkBinaryThinningImageFilter3D.h itkBinaryThinningImageFilter3D.txx)

ADD_LIBRARY(helpers helpers.cpp helpers.h)
ADD_LIBRARY(ftkTrackFeatures ftkTrackFeatures.h ftkTrackFeatures.cpp)

ADD_EXECUTABLE(unmix unmix.cpp ${SRCS})
ADD_EXECUTABLE(segmentation segmentation.cpp ${SRCS})
ADD_EXECUTABLE(tracking tracking.cpp ${SRCS})
ADD_EXECUTABLE(vessel_thinning vessel_thinning.cpp ${SRCS})
ADD_EXECUTABLE(summary summary.cpp ${SRCS})


TARGET_LINK_LIBRARIES(unmix ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers)
TARGET_LINK_LIBRARIES(segmentation ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers)
TARGET_LINK_LIBRARIES(tracking ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers)
TARGET_LINK_LIBRARIES(vessel_thinning ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers)
TARGET_LINK_LIBRARIES(summary ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers ftkTrackFeatures)

INSTALL(TARGETS unmix RUNTIME DESTINATION bin)
INSTALL(TARGETS segmentation RUNTIME DESTINATION bin)
INSTALL(TARGETS tracking RUNTIME DESTINATION bin)
INSTALL(TARGETS vessel_thinning RUNTIME DESTINATION bin)
INSTALL(TARGETS summary RUNTIME DESTINATION bin)

