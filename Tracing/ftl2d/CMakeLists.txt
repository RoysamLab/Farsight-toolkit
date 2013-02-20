# CMakeList.txt for ftl2d Farsight
# Author: Amit Mukherjee
PROJECT ( ftl2d_SETracing )

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
             INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
             MESSAGE(FATAL_ERROR
             "Cannot build SETracing without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(LibXml2)
IF(LIBXML2_FOUND)
             INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
ELSE(LIBXML2_FOUND)
             MESSAGE(FATAL_ERROR
             "Cannot build SETracing without LIBXML2.  Please set LIBXML2_DIR.")
ENDIF(LIBXML2_FOUND)



SET ( ftl2d_SETracing_sources
      ftl2d_SEtracing2D.cxx
      ftl2d_TraceConfig.cxx
      ftl2d_SeedContainer2D.cxx
      ftl2d_Segment2D.cxx
      ftl2d_SuperEllipsoid2D.cxx
      ftl2d_Tracer.cxx
      ftl2d_Vessel.cxx
      ftl2d_SeedPoint2D.cxx
)

SET ( ftl2d_SETracing_headers 
      ftl2d_TraceConfig.h
      ftl2d_SeedContainer2D.h
      ftl2d_SeedPoint2D.h
      ftl2d_Tracer.h
      ftl2d_Segment2D.h
      ftl2d_SuperEllipsoid2D.h
      ftl2d_Vessel.h
      ftl2d_BranchPoint.h
)     


ADD_EXECUTABLE ( ftl2d_SETracing ${ftl2d_SETracing_sources} ${ftl2d_SETracing_headers} )
TARGET_LINK_LIBRARIES( ftl2d_SETracing 
  ${ITK_LIBRARIES} 
  ${LIBXML2_LIBRARIES}
)
