################################################################################
################################################################################


SET( FTKNUCLEARSEGMENTATIONNIC_SRCS

main.cpp
#1__1_1_BinarizeMixPoisson/graph.cpp
#1__1_1_BinarizeMixPoisson/maxflow.cpp
#ftkNuclearSegmentationNic.cpp
 
)

################################################################################
################################################################################


SET( FTKNUCLEARSEGMENTATIONNIC_HDRS

ftkGlobalIncludes.h
frkGlobalStructs.h

ftkNuclearSegmentationNic.h
ftkNuclearSegmentationNic.hxx

BinarizeMixPoisson/BinarizeMixPoisson.h
BinarizeMixPoisson/graph.h

SeedDetection/SeedDetectionLoG.h
SeedDetection/itklaplacianrecursivegaussianimagefilternew.h

MaxClustering/MaxClustering.h
 
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
