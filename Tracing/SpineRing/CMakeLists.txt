SET( Libraries
	${ITK_LIBRARIES}
	TinyXML
)



SET ( SpineSegmentation_sources
      RingSE.cxx
      SpineRing.cxx
      TraceSegNode.cxx
      ImageProc.cxx
      SpinePaths.cxx
      FMPath.cxx
      SpeedFunc.cxx

      itkSingleImageCostFunction.h
      itkArrivalFunctionToPathFilter.h
      itkSpeedFunctionToPathFilter.h
      SpineRing.h
      SpineConsts.h
      SpineUtils.h
      TraceSegNode.h
      CommonTypeDefs.h
      SpineTypedefs.h
      ImageProc.h
)      


INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}
  ${SOURCE_PATH}
 )

ADD_EXECUTABLE ( SpineSegmentation ${SpineSegmentation_sources} SpineMain.cxx )
TARGET_LINK_LIBRARIES( SpineSegmentation ${Libraries} )
INSTALL(TARGETS SpineSegmentation DESTINATION ${INSTALL_BIN_DIR} )

