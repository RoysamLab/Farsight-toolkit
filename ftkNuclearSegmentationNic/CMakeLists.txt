ADD_SUBDIRECTORY(

#1__1_1_BinarizeMixPoisson
#1__1_2_BinarizeMixGaussian
#1__1_3_BinarizeOtsu
#1__1_4_OtherMethod
#1__2_1_SeedDetection
#1__3_1_MaxClustering
#1__4_1_AlphaExpansion
#2__1_1_RadVoting
#2__2_1_CellShapeIni
#2__3_1_GVFforce
#2__3_2_AnotherForce
#2__4_1_LevelSet
#3__1_1_LevlSet

#)

################################################################################
################################################################################


SET( FTKNUCLEARSEGMENTATIONNIC_SRCS

main.cpp
1__1_1_BinarizeMixPoisson/graph.cpp
1__1_1_BinarizeMixPoisson/maxflow.cpp
#ftkNuclearSegmentationNic.cpp
 
)

################################################################################
################################################################################


SET( FTKNUCLEARSEGMENTATIONNIC_HDRS

ftkGlobalIncludes.h
frkGlobalStructs.h

ftkNuclearSegmentationNic.h
ftkNuclearSegmentationNic.hxx

1__1_1_BinarizeMixPoisson/BinarizeMixPoisson.h
1__1_1_BinarizeMixPoisson/graph.h

1__2_1_SeedDetection/SeedDetectionLoG.h
1__2_1_SeedDetection/itklaplacianrecursivegaussianimagefilternew.h

1__3_1_MaxClustering/MaxClustering.h
 
)

################################################################################
################################################################################


ADD_LIBRARY( libftkNuclearSegmentationNic ${FTKNUCLEARSEGMENTATIONNIC_HDRS} ${FTKNUCLEARSEGMENTATIONNIC_SRCS})

ADD_EXECUTABLE( ftkNuclearSegmentationNic
		main.cpp
)

################################################################################
################################################################################


TARGET_LINK_LIBRARIES ( ftkNuclearSegmentationNic ${VTK_LIBRARIES} ${ITK_LIBRARIES} ftkCommon libftkNuclearSegmentationNic )
