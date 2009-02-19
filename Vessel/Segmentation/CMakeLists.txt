ADD_EXECUTABLE(vessel_segmentation main.cpp find_median.cpp)
TARGET_LINK_LIBRARIES(vessel_segmentation ITKCommon ITKIO )

INSTALL(TARGETS vessel_segmentation DESTIONATION bin)
