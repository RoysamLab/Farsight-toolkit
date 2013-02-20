PROJECT(Curvelets)

INCLUDE(FindFFTW215.cmake)
ADD_LIBRARY(CurveletClass
		Curvelet.h	Curvelet.cpp
		)
		
		QT4_WRAP_CPP( GUI_MOC_SRCS CurveletGUI.h)
ADD_LIBRARY(CurveletClassGUI
		CurveletGUI.h	CurveletGUI.cpp
  ${GUI_MOC_SRCS} 
		)
IF(FFTW_FOUND)
	INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
	ADD_LIBRARY(fdct_wrapping STATIC fdct_wrapping.cpp ifdct_wrapping.cpp fdct_wrapping_param.cpp)
	TARGET_LINK_LIBRARIES(fdct_wrapping ${FFTW_LIBRARIES})
	ADD_EXECUTABLE(curvelets curvelet_preprocessing.cpp)
	TARGET_LINK_LIBRARIES(curvelets CurveletClass fdct_wrapping        ${ITK_LIBRARIES} )
	INSTALL(TARGETS curvelets RUNTIME DESTINATION ${INSTALL_BIN_DIR})
	
	ADD_EXECUTABLE(CurveletGUI CurveletGUIDriver.cpp)
	TARGET_LINK_LIBRARIES(CurveletGUI CurveletClassGUI CurveletClass fdct_wrapping    ${QT_LIBRARIES}      ${ITK_LIBRARIES} )
	INSTALL(TARGETS CurveletGUI RUNTIME DESTINATION ${INSTALL_BIN_DIR})
	
	INSTALL(TARGETS fdct_wrapping
	  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
	  LIBRARY DESTINATION ${INSTALL_LIB_DIR} )
ENDIF(FFTW_FOUND)
