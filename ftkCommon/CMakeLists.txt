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
    SET ( FTKCOMMON_SRCS ${FTKCOMMON_SRCS} ftkLaplacianOfGaussian3DCPU.cpp )
    SET ( FTKCOMMON_HDRS ${FTKCOMMON_HDRS} ftkLaplacianOfGaussian3DCPU.h)
endif (BUILD_CUDA)


# GnT -> Graphs and Tree Classes
# Contains RAG and Merge Tree Codes

SET(GnT_SRCS
ftkgnt.cpp
)

SET(GnT_HDRS
ftkgnt.h
)


FIND_PACKAGE(Boost)

LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS})
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS})

ADD_LIBRARY( ftkCommon
            ${FTKCOMMON_SRCS} 
            ${FTKCOMMON_HDRS} 
            ${GnT_HDRS}
            ${GnT_SRCS}
)


TARGET_LINK_LIBRARIES( ftkCommon TinyXML ftkImage
                       vtksys vtkCommon vtkFiltering vtkIO
                       ITKCommon ITKBasicFilters ITKIO ITKStatistics )

ADD_EXECUTABLE( block main.cpp )
TARGET_LINK_LIBRARIES( block ftkCommon )

INSTALL( TARGETS block RUNTIME DESTINATION bin )

