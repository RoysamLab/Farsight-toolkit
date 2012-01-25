OPTION( BUILD_MDL "Build MDL Tracing" OFF )

SUBDIRS( Focus )
SUBDIRS(MontagingToTraceEditorXMLConverter)

IF( BUILD_MDL )
    #SUBDIRS( MDL )
    SUBDIRS( MDL2 )
ENDIF()

OPTION( BUILD_MDL_XL "Build MDL_XL" OFF )
IF( BUILD_MDL_XL )
    SUBDIRS( MDL_XL )
ENDIF()

OPTION( BUILD_ftl2d "Build FTL2D Tracing" OFF )
IF( BUILD_ftl2d )
  SUBDIRS( ftl2d )
ENDIF()

OPTION( BUILD_TraceEdit "Build Trace Editor" ON )
IF( BUILD_TraceEdit )
  SUBDIRS( TraceEdit )
ENDIF()

OPTION( BUILD_ftl3d "Build FTL3D Tracing" OFF )
IF( BUILD_ftl3d )
    SUBDIRS( ftl3d )
ENDIF()

OPTION( BUILD_RPITrace3D "Build RPI Trace3D" OFF )
IF( BUILD_RPITrace3D )
  SUBDIRS( RPITrace3D )
ENDIF()

OPTION( BUILD_SpineRing "Build Spine Ring" OFF )
IF( BUILD_SpineRing )
  SUBDIRS( SpineRing )
ENDIF()

OPTION( BUILD_Preprocessing "Build the fast tensor voting library and executable" OFF )
IF( BUILD_Preprocessing )
  SUBDIRS( Preprocessing/SHtensor_cpp/2D )
ENDIF()

OPTION( BUILD_Rendering "Build 5-Label Result Renderer" ON )
IF( BUILD_Rendering )
  SUBDIRS( Rendering )
ENDIF()

OPTION( BUILD_ExtractSomas "Build Soma Extraction executables" ${ITK_USE_REVIEW} )
IF( BUILD_ExtractSomas )
  SUBDIRS( ExtractSomas )
ENDIF()

OPTION( BUILD_MultipleNeuronTracer "Build Multiple Neuron Tracer" ON )
IF( BUILD_MultipleNeuronTracer )
  SUBDIRS( MultipleNeuronTracer )
ENDIF()

OPTION( BUILD_Microglia "Build Microglia" OFF )
IF( BUILD_Microglia )
    SUBDIRS( Microglia )
ENDIF()


OPTION( BUILD_TracingSystem "Build The Tracing System" OFF )
IF(BUILD_TracingSystem)
  SUBDIRS( TheTracingSystem )
  SET(CPACK_PACKAGE_EXECUTABLES "TheTracingSystem" "The Tracing System" ${CPACK_PACKAGE_EXECUTABLES})
ENDIF()


OPTION( BUILD_MicrogliaRegionTracer "Build Microglia Region Tracer" ON )
IF( BUILD_MicrogliaRegionTracer )
  SUBDIRS( MicrogliaRegionTracer )
ENDIF()



