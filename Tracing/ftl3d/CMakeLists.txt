 # CMakeList.txt for  Trace3D
# Author: Amit Mukherjee
PROJECT(TraceSE3D)
cmake_minimum_required( VERSION 2.6)
if(COMMAND cmake_policy)
   cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET( Sources Trace3dDriver.cxx TraceMain.cxx TraceMain.h)

ADD_LIBRARY(TraceSE3DLib
  TraceConfig.cxx       TraceConfig.h
  TVessel.cxx           TVessel.h
  SegInit.cxx           SegInit.h
  SegFit.cxx            SegFit.h
  SeedPoint3D.cxx       SeedPoint3D.h
  SeedContainer3D.cxx   SeedContainer3D.h       
  Trace.cxx             Trace.h
  Seed2Seg.cxx          Seed2Seg.h
  NodeContainer3D.cxx   NodeContainer3D.h
  TraceContainer3D.cxx  TraceContainer3D.h
  TraceNode.h
)      
TARGET_LINK_LIBRARIES( TraceSE3DLib  TinyXML      ${ITK_LIBRARIES} )

SET ( Libraries
      ${ITK_LIBRARIES}
      TraceSE3DLib
)
QT4_WRAP_CPP( GUI_MOC_SRCS TraceMain.h)

ADD_EXECUTABLE( ${PROJECT_NAME} ${EXE_TYPE} ${Sources} ${GUI_MOC_SRCS} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${QT_LIBRARIES} ${Libraries} ftkCommon )
SET(CPACK_PACKAGE_EXECUTABLES "${PROJECT_NAME};Trace FL3D")

ADD_EXECUTABLE ( RenderTraces Render.cxx TraceNode.h ${XML_sources})
TARGET_LINK_LIBRARIES( RenderTraces
                       ${Libraries}
                       ${ITK_LIBRARIES}
                       tiff
		       ${VTK_LIBRARIES}
)

INSTALL( TARGETS RenderTraces DESTINATION ${INSTALL_BIN_DIR} )
