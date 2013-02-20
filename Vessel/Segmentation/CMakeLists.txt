PROJECT (vessel_segmentation)
ADD_EXECUTABLE(vessel_segmentation main.cpp )

SET(vessel_segmentation_sources
	find_median.cpp
    itkHessianToObjectnessMeasureImageFilter.txx	
	itkMultiScaleHessianBasedMeasureImageFilter.txx
)

SET(vessel_segmentation_headers
	itkHessianToObjectnessMeasureImageFilter.h
	itkMultiScaleHessianBasedMeasureImageFilter.h
)

TARGET_LINK_LIBRARIES(vessel_segmentation ${ITK_LIBRARIES} )

TARGET_LINK_LIBRARIES(vessel_segmentation ${ITK_LIBRARIES} ${VXL_LIBRARIES} rgtl rtvl)

INSTALL(TARGETS vessel_segmentation DESTINATION ${INSTALL_BIN_DIR} )
