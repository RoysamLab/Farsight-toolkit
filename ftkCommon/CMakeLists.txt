SET( FTKCOMMON_SRCS
  ftkUtils.cpp
  ftkPreprocess.cpp
  ftkMultipleImageHandler.cpp
  ftkParameters.cpp
  ftkProjectManager.cxx
ftkLaplacianOfGaussian3D.cpp
)

SET( FTKCOMMON_HDRS
  ftkUtils.h
  ftkPreprocess.h
  ftkMultipleImageHandler.h
  ftkParameters.h
  ftkProjectManager.h
  ftkLaplacianOfGaussian3D.h
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
                       vtksys vtkCommon vtkFiltering vtkIO
                       ${ITK_LIBRARIES} )

ADD_EXECUTABLE( block main.cpp )
TARGET_LINK_LIBRARIES( block ftkCommon )

INSTALL( TARGETS block RUNTIME DESTINATION bin )

#16 bit to 8 bit TIF convertor
ADD_EXECUTABLE ( TIF16bitTo8bit TIF16bitTo8Bit.cpp )
TARGET_LINK_LIBRARIES( TIF16bitTo8bit      ${ITK_LIBRARIES} )
INSTALL( TARGETS TIF16bitTo8bit RUNTIME DESTINATION bin )

#ADD_EXECUTABLE(ftkLoG3D_tester ftkLoG3D_tester.cpp)
#TARGET_LINK_LIBRARIES(ftkLoG3D_tester ${ITK_LIBRARIES} ftkCommon)
#INSTALL(TARGETS ftkLoG3D_tester RUNTIME DESTINATION bin)
