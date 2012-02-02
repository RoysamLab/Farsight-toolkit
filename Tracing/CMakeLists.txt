OPTION( BUILD_MDL "Build MDL Tracing" OFF )

ADD_SUBDIRECTORY( Focus )
ADD_SUBDIRECTORY(MontagingToTraceEditorXMLConverter)

IF( BUILD_MDL )
    #ADD_SUBDIRECTORY( MDL )
    ADD_SUBDIRECTORY( MDL2 )
ENDIF()

OPTION( BUILD_MDL_XL "Build MDL_XL" OFF )
IF( BUILD_MDL_XL )
    ADD_SUBDIRECTORY( MDL_XL )
ENDIF()

OPTION( BUILD_ftl2d "Build FTL2D Tracing" OFF )
IF( BUILD_ftl2d )
  ADD_SUBDIRECTORY( ftl2d )
ENDIF()

OPTION( BUILD_TraceEdit "Build Trace Editor" ON )
IF( BUILD_TraceEdit )
  ADD_SUBDIRECTORY( TraceEdit )
ENDIF()

OPTION( BUILD_ftl3d "Build FTL3D Tracing" OFF )
IF( BUILD_ftl3d )
    ADD_SUBDIRECTORY( ftl3d )
ENDIF()

OPTION( BUILD_RPITrace3D "Build RPI Trace3D" OFF )
IF( BUILD_RPITrace3D )
  ADD_SUBDIRECTORY( RPITrace3D )
ENDIF()

OPTION( BUILD_SpineRing "Build Spine Ring" OFF )
IF( BUILD_SpineRing )
  ADD_SUBDIRECTORY( SpineRing )
ENDIF()

OPTION( BUILD_Preprocessing "Build the fast tensor voting library and executable" OFF )
IF( BUILD_Preprocessing )
  ADD_SUBDIRECTORY( Preprocessing/SHtensor_cpp/2D )
ENDIF()

OPTION( BUILD_Rendering "Build 5-Label Result Renderer" ON )
IF( BUILD_Rendering )
  ADD_SUBDIRECTORY( Rendering )
ENDIF()

OPTION( BUILD_ExtractSomas "Build Soma Extraction executables" ${ITK_USE_REVIEW} )
IF( BUILD_ExtractSomas )
  ADD_SUBDIRECTORY( ExtractSomas )
ENDIF()

OPTION( BUILD_MultipleNeuronTracer "Build Multiple Neuron Tracer" ON )
IF( BUILD_MultipleNeuronTracer )
  ADD_SUBDIRECTORY( MultipleNeuronTracer )
ENDIF()

OPTION( BUILD_Microglia "Build Microglia" OFF )
IF( BUILD_Microglia )
    ADD_SUBDIRECTORY( Microglia )
ENDIF()


OPTION( BUILD_TracingSystem "Build The Tracing System" OFF )
IF(BUILD_TracingSystem)
  ADD_SUBDIRECTORY( TheTracingSystem )
  SET(CPACK_PACKAGE_EXECUTABLES "TheTracingSystem" "The Tracing System" ${CPACK_PACKAGE_EXECUTABLES})
ENDIF()


OPTION( BUILD_MicrogliaRegionTracer "Build Microglia Region Tracer" OFF )
IF( BUILD_MicrogliaRegionTracer )
  ADD_SUBDIRECTORY( MicrogliaRegionTracer )
ENDIF()



