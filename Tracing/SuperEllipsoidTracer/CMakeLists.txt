# CMakeList.txt for  Trace3D
# Author: Amit Mukherjee

cmake_minimum_required( VERSION 2.6)
if(COMMAND cmake_policy)
   cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


PROJECT ( Trace3D )

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
             INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
             MESSAGE(FATAL_ERROR
             "Cannot build SETracing without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
             INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
             MESSAGE(FATAL_ERROR "Cannot build without VTK.  Please set VTK_DIR.")
ENDIF(VTK_FOUND)

SET (XML_sources
	tinystr.cpp
	tinyxml.cpp
	tinyxmlerror.cpp
	tinyxmlparser.cpp
	tinystr.h
	tinyxml.h
)	



SET ( Trace3D_sources
      TraceConfig.cxx
      TVessel.cxx
      SegInit.cxx
      SegFit.cxx
      SeedPoint3D.cxx
      SeedContainer3D.cxx	      
      Trace.cxx
      SegInit.cxx
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

#ADD_LIBRARY ( TinyXML ${XML_sources})

ADD_EXECUTABLE ( Trace3D ${Trace3D_sources} ${Trace3D_headers}  TraceMain.cxx )
#TARGET_LINK_LIBRARIES( Trace3D ITKCommon ITKIO ITKNumerics TinyXML)
TARGET_LINK_LIBRARIES( Trace3D ${ITK_LIBRARIES} TinyXML) 

ADD_EXECUTABLE ( RenderTraces Render.cxx TraceNode.h )
TARGET_LINK_LIBRARIES( RenderTraces vtkRendering vtkIO 
vtkCommon  vtksys vtkFiltering  vtkVolumeRendering vtkWidgets TinyXML ${ITK_LIBRARIES})

ADD_EXECUTABLE ( ConvertToSWC ConvertToSWC.cxx TraceNode.h )
TARGET_LINK_LIBRARIES( ConvertToSWC TinyXML ${ITK_LIBRARIES})

ADD_EXECUTABLE ( TraceOverlay TraceOverlay.cxx TraceNode.h )
TARGET_LINK_LIBRARIES( TraceOverlay TinyXML ${ITK_LIBRARIES})


ADD_EXECUTABLE( Trace3D_package ${Trace3D_sources} ${Trace3D_headers} TraceMain.cxx Render.cxx ConvertToSWC.cxx TraceOverlay.cxx TraceNode.h )
TARGET_LINK_LIBRARIES( Trace3D_package ${ITK_LIBRARIES} TinyXML vtkRendering vtkIO vtkCommon  vtksys vtkFiltering  vtkVolumeRendering vtkWidgets )
install(TARGETS Trace3D_package DESTINATION ${INSTALL_BIN_DIR} )