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




TARGET_LINK_LIBRARIES(vessel_segmentation        ${ITK_LIBRARIES} )
include_directories(${VXL_INCLUDE_DIR})

FIND_PACKAGE(VXL)
IF(VXL_FOUND)
    INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
ENDIF(VXL_FOUND) 

IF(VXL_FOUND)
    #INCLUDE(${VXL_RPL_INCLUDE_DIR}/UseVXL.cmake)
ENDIF(VXL_FOUND)



IF(VXL_FOUND)
	#TARGET_LINK_LIBRARIES(vessel_segmentation ${ITK_LIBRARIES} rgtl rtvl)
	TARGET_LINK_LIBRARIES(vessel_segmentation ${ITK_LIBRARIES} ${VXL_LIBRARIES} rgtl rtvl)
ELSE(VXL_FOUND)
	 MESSAGE(FATAL_ERROR "Cannot build without  VXL.")
ENDIF(VXL_FOUND)


INSTALL(TARGETS vessel_segmentation DESTINATION ${INSTALL_BIN_DIR} )
