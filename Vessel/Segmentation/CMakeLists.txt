ADD_EXECUTABLE(vessel_segmentation main.cpp )
TARGET_LINK_LIBRARIES(vessel_segmentation ITKCommon ITKIO )

INSTALL(TARGETS vessel_segmentation DESTINATION ${INSTALL_BIN_DIR} )
