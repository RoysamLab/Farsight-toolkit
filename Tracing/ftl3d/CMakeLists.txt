 # CMakeList.txt for  Trace3D
# Author: Amit Mukherjee

cmake_minimum_required( VERSION 2.6)
if(COMMAND cmake_policy)
   cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

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

SET ( Trace3D_libraries
      ITKCommon
      ITKIO
      ITKNumerics
      TinyXML   
)


ADD_EXECUTABLE ( Trace3D ${Trace3D_sources} ${Trace3D_headers}  TraceMain.cxx )
TARGET_LINK_LIBRARIES( Trace3D  ${Trace3D_libraries})

ADD_EXECUTABLE ( RenderTraces Render.cxx TraceNode.h ${XML_sources})
TARGET_LINK_LIBRARIES( RenderTraces
                       ${Trace3D_libraries}
                       tiff
                       ITKBasicFilters
                       vtkRendering
                       vtkIO 
                       vtkCommon
                       vtksys
                       vtkFiltering
 		               vtkVolumeRendering 
					   vtkWidgets 
					   TinyXML
                        )

ADD_EXECUTABLE ( ConvertToSWC ConvertToSWC.cxx TraceNode.h )
TARGET_LINK_LIBRARIES( ConvertToSWC ITKCommon ITKBasicFilters ITKIO ITKNumerics TinyXML )

INSTALL( TARGETS Trace3D DESTINATION bin )
INSTALL( TARGETS RenderTraces DESTINATION bin )
INSTALL( TARGETS ConvertToSWC DESTINATION bin )

