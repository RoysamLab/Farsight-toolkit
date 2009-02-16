SET( NSRCS
	yousef_seg.h
	yousef_seg.cpp
)


################################################################################
#      Here we add all sub-directories that are involved with this program       #
################################################################################

SUBDIRS( cell_binarization )
SUBDIRS( seed_detection )
SUBDIRS( seed_detection_3D )
SUBDIRS( clustering_3D )
SUBDIRS( graphColLearn_3D )
SUBDIRS( alpha_expansion_3d )
SUBDIRS( cellLabelToBorders )
SUBDIRS( EM_GMM )


ADD_LIBRARY( Yousef_Nucleus_Seg ${NSRCS} )

TARGET_LINK_LIBRARIES ( Yousef_Nucleus_Seg 
						cell_binarization
						seedsDetection_2D
						seedsDetection_3D
						clustering_3D
						graphColLearn_3D
						alpha_expansion_3d
						cellLabelToBorders
						EM_GMM
						ITKBasicFilters ITKCommon ITKIO)
