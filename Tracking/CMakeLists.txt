PROJECT(Tracker)

cmake_minimum_required(VERSION 2.0)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)



SET(SRCS helpers.h itkBinaryThinningImageFilter3D.h ftkTrackFeatures.h itkBinaryThinningImageFilter3D.h itkBinaryThinningImageFilter3D.txx)

ADD_LIBRARY(helpers helpers.cpp helpers.h)
ADD_LIBRARY(ftkTrackFeatures ftkTrackFeatures.h ftkTrackFeatures.cpp)

ADD_EXECUTABLE(unmix unmix.cpp helpers.h)
ADD_EXECUTABLE(segmentation segmentation.cpp helpers.h)
ADD_EXECUTABLE(tracking tracking.cpp helpers.h ftkTrackFeatures.h)
ADD_EXECUTABLE(vessel_thinning vessel_thinning.cpp ${SRCS})
ADD_EXECUTABLE(summary summary.cpp helpers.h ftkTrackFeatures.h)
ADD_EXECUTABLE(background_subtraction background_subtraction.cpp helpers.h)
ADD_EXECUTABLE(extract_slices extract_slices.cpp helpers.h)
ADD_EXECUTABLE(anisodiff anisodiff.cpp helpers.h)

TARGET_LINK_LIBRARIES(unmix ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers)
TARGET_LINK_LIBRARIES(segmentation ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers)
TARGET_LINK_LIBRARIES(tracking ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers)
TARGET_LINK_LIBRARIES(vessel_thinning ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers)
TARGET_LINK_LIBRARIES(summary ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers ftkTrackFeatures)
TARGET_LINK_LIBRARIES(background_subtraction ITKCommon ITKIO itkvnl itkvnl_algo Yousef_Nucleus_Seg ITKStatistics helpers)
TARGET_LINK_LIBRARIES(extract_slices ITKCommon ITKIO itkvnl itkvnl_algo ITKStatistics helpers)
TARGET_LINK_LIBRARIES(anisodiff ITKCommon ITKIO itkvnl itkvnl_algo ITKStatistics helpers)

INSTALL(TARGETS unmix RUNTIME DESTINATION bin)
INSTALL(TARGETS segmentation RUNTIME DESTINATION bin)
INSTALL(TARGETS tracking RUNTIME DESTINATION bin)
INSTALL(TARGETS vessel_thinning RUNTIME DESTINATION bin)
INSTALL(TARGETS summary RUNTIME DESTINATION bin)
INSTALL(TARGETS background_subtraction RUNTIME DESTINATION bin)
INSTALL(TARGETS extract_slices RUNTIME DESTINATION bin)
INSTALL(TARGETS anisodiff RUNTIME DESTINATION bin)

