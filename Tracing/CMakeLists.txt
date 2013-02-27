option( BUILD_MDL "Build MDL Tracing" OFF )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

add_subdirectory( Focus )
add_subdirectory( MontagingToTraceEditorXMLConverter )

if( BUILD_MDL )
  #add_subdirectory( MDL )
  add_subdirectory( MDL2 )
endif()

option( BUILD_MDL_XL "Build MDL_XL" OFF )
if( BUILD_MDL_XL )
  add_subdirectory( MDL_XL )
endif()

option( BUILD_ftl2d "Build FTL2D Tracing" OFF )
if( BUILD_ftl2d )
  add_subdirectory( ftl2d )
endif()

option( BUILD_TraceEdit "Build Trace Editor" ON )
if( BUILD_TraceEdit )
  add_subdirectory( TraceEdit )
endif()

option( BUILD_ftl3d "Build FTL3D Tracing" OFF )
if( BUILD_ftl3d )
    add_subdirectory( ftl3d )
endif()

option( BUILD_RPITrace3D "Build RPI Trace3D" OFF )
if( BUILD_RPITrace3D )
  add_subdirectory( RPITrace3D )
endif()

option( BUILD_SpineRing "Build Spine Ring" OFF )
if( BUILD_SpineRing )
  add_subdirectory( SpineRing )
endif()

option( BUILD_Preprocessing "Build the fast tensor voting library and executable" OFF )
if( BUILD_Preprocessing )
  add_subdirectory( Preprocessing/SHtensor_cpp/2D )
endif()

option( BUILD_Rendering "Build 5-Label Result Renderer" OFF )
if( BUILD_Rendering )
  add_subdirectory( Rendering )
endif()

option( BUILD_MultipleNeuronTracer "Build Multiple Neuron Tracer" OFF )
if( BUILD_MultipleNeuronTracer )
  add_subdirectory( MultipleNeuronTracer )
endif()

option( BUILD_Microglia "Build Microglia Tracing Module" OFF )
if( BUILD_Microglia )
    add_subdirectory( Microglia )
endif()

option( BUILD_TracingSystem "Build The Tracing System" OFF )
if(BUILD_TracingSystem)
  add_subdirectory( TheTracingSystem )
  set(CPACK_PACKAGE_EXECUTABLES "TheTracingSystem" "The Tracing System" ${CPACK_PACKAGE_EXECUTABLES})
endif()

option( BUILD_MicrogliaRegionTracer "Build Microglia Region Tracer" OFF )
if( BUILD_MicrogliaRegionTracer )
  add_subdirectory( MicrogliaRegionTracer )
endif()

option( BUILD_AstroTracer "Build Astrocyte Tracer" OFF )
if( BUILD_AstroTracer )
  add_subdirectory( AstroTracer )
endif()

option( BUILD_SuperEllipsoidTracer "Build SuperEllipsoid Tracer" OFF )
if( BUILD_SuperEllipsoidTracer )
  add_subdirectory( SuperEllipsoidTracer )
endif()

option( BUILD_VESSEL_TRACING "Build Ball Tracing for Vessels" OFF )
if( BUILD_VESSEL_TRACING )
  set(CPACK_PACKAGE_EXECUTABLES "ftkVesselTracer" "Vessel Tracer" ${CPACK_PACKAGE_EXECUTABLES})
  add_subdirectory( ftkVesselTracer )
  add_definitions( -DUSE_BALL_TRACER )
  include_directories( ${PROJECT_SOURCE_DIR}/Tracing/ftkVesselTracer )
endif()