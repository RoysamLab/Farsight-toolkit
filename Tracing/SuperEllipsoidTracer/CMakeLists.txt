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

ADD_EXECUTABLE ( Trace3D ${Trace3D_sources} ${Trace3D_headers}  TraceMain.cxx )
TARGET_LINK_LIBRARIES( Trace3D ${ITK_LIBRARIES} TinyXML) 

ADD_EXECUTABLE ( RenderTracesSuperEllipsoid Render.cxx TraceNode.h )
TARGET_LINK_LIBRARIES( RenderTracesSuperEllipsoid ${VTK_LIBRARIES} TinyXML ${ITK_LIBRARIES})

ADD_EXECUTABLE ( ConvertToSWC ConvertToSWC.cxx TraceNode.h )
TARGET_LINK_LIBRARIES( ConvertToSWC TinyXML ${ITK_LIBRARIES})

ADD_EXECUTABLE ( TraceOverlay TraceOverlay.cxx TraceNode.h )
TARGET_LINK_LIBRARIES( TraceOverlay TinyXML ${ITK_LIBRARIES})

#ADD_EXECUTABLE( Trace3D_package ${Trace3D_sources} ${Trace3D_headers} TraceMain.cxx Render.cxx ConvertToSWC.cxx TraceOverlay.cxx TraceNode.h )
#TARGET_LINK_LIBRARIES( Trace3D_package ${ITK_LIBRARIES} TinyXML ${VTK_LIBRARIES} )
