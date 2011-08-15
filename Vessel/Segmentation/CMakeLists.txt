ADD_EXECUTABLE(vessel_segmentation main.cpp )
TARGET_LINK_LIBRARIES(vessel_segmentation        ${ITK_LIBRARIES} )

INSTALL(TARGETS vessel_segmentation DESTINATION ${INSTALL_BIN_DIR} )
