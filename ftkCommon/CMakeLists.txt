SET( FTKCOMMON_SRCS
  ftkUtils.cpp
  ftkPreprocess.cpp
  ftkMultipleImageHandler.cpp
  ftkParameters.cpp
  ftkProjectManager.cxx
)

SET( FTKCOMMON_HDRS
  ftkUtils.h
  ftkPreprocess.h
  ftkMultipleImageHandler.h
  ftkParameters.h
  ftkProjectManager.h
)

if (BUILD_CUDA)
    SET ( FTKCOMMON_SRCS ${FTKCOMMON_SRCS} Convolution_CUDA.cu )
    SET ( FTKCOMMON_HDRS ${FTKCOMMON_HDRS} Convolution_CUDA.cuh )
endif (BUILD_CUDA)


# GnT -> Graphs and Tree Classes
# Contains RAG and Merge Tree Codes

SET(GnT_SRCS
ftkgnt.cpp
)

SET(GnT_HDRS
ftkgnt.h
)


if (BUILD_CUDA)
	CUDA_ADD_LIBRARY( ftkCommon
            ${FTKCOMMON_SRCS} 
            ${FTKCOMMON_HDRS} 
            ${GnT_HDRS}
            ${GnT_SRCS}
	)
else (BUILD_CUDA)
	ADD_LIBRARY( ftkCommon
            ${FTKCOMMON_SRCS} 
            ${FTKCOMMON_HDRS} 
            ${GnT_HDRS}
            ${GnT_SRCS}
	)
endif (BUILD_CUDA)


TARGET_LINK_LIBRARIES( ftkCommon TinyXML ftkImage QHULL
                       ${VTK_LIBRARIES}
                       ${ITK_LIBRARIES} )

ADD_EXECUTABLE( block main.cpp )
TARGET_LINK_LIBRARIES( block ftkCommon )

INSTALL( TARGETS block RUNTIME DESTINATION ${INSTALL_BIN_DIR} )

#16 bit to 8 bit TIF convertor
ADD_EXECUTABLE ( TIF16bitTo8bit TIF16bitTo8Bit.cpp )
TARGET_LINK_LIBRARIES( TIF16bitTo8bit      ${ITK_LIBRARIES} )
INSTALL( TARGETS TIF16bitTo8bit RUNTIME DESTINATION ${INSTALL_BIN_DIR} )

ADD_EXECUTABLE(itkRolloverOpenMPTest itkRolloverOpenMPTest.cpp)
TARGET_LINK_LIBRARIES(itkRolloverOpenMPTest ${ITK_LIBRARIES})

#Removing itkRolloverOpenMPTest because it takes way too long to run on dashboard machines
#ADD_TEST(itkRolloverOpenMPTest
#  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkRolloverOpenMPTest
#  ${FARSIGHT_DATA_ROOT}/TraceData/10_spine2.tif
#  )
#SET_TESTS_PROPERTIES(itkRolloverOpenMPTest PROPERTIES
#  TIMEOUT 1800) #30 minutes

ADD_EXECUTABLE(ftkVS2010CrashTest ftkVS2010CrashTest.cpp)
TARGET_LINK_LIBRARIES(ftkVS2010CrashTest ${ITK_LIBRARIES} ${QT_LIBRARIES})
ADD_TEST(ftkVS2010Crashtest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ftkVS2010CrashTest
  ${FARSIGHT_DATA_ROOT}/TraceData/10_spine2.tif
  )
