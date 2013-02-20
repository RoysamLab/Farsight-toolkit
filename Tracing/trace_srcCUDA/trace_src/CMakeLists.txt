# CMakeList.txt for  Trace3D
# Author: Amit Mukherjee
cmake_minimum_required( VERSION 2.4)

PROJECT ( Trace3D )

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/cuda" ${CMAKE_MODULE_PATH})

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
             INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
             MESSAGE(FATAL_ERROR
             "Cannot build SETracing without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(CUDA REQUIRED)


SET ( Trace3D_sources
      TraceConfig.cxx
      TVessel.cxx
      SegFit.cxx
      SeedPoint3D.cxx
      SeedContainer3D.cxx	      
      Trace.cxx
      Seed2Seg.cxx
      NodeContainer3D.cxx
      TraceContainer3D.cxx
)

SET ( Trace3D_headers 
      TraceConfig.h
      TVessel.h
      SeedContainer3D.h
      SegInit.h
      SegFit.h
      SeedPoint3D.h
      Trace.h
      Seed2Seg.h
      NodeContainer3D.h
      TraceContainer3D.h
)      

CUDA_ADD_EXECUTABLE ( Trace3D SegInit.cu ${Trace3D_sources} ${Trace3D_headers} TraceMain.cxx )
TARGET_LINK_LIBRARIES( Trace3D        ${ITK_LIBRARIES} )

