SET( NSRCS
	yousef_seg.h
	yousef_seg.cpp
)

################################################################################
#      Here we add all sub-directories that are involved with this program       #
################################################################################

#SUBDIRS( cell_binarization )
#SUBDIRS( seed_detection )
#SUBDIRS( seed_detection_3D )
#SUBDIRS( clustering_3D )
#SUBDIRS( graphColLearn_3D )
#SUBDIRS( alpha_expansion_3d )
#SUBDIRS( cellLabelToBorders )
#SUBDIRS( EM_GMM )

SET( sd3d_SRCS
	seed_detection_3D/seedsDetection_3D.cxx
	seed_detection_3D/seedsdetection_3D.h
	seed_detection_3D/itklaplacianrecursivegaussianimagefilternew.txx
	seed_detection_3D/itklaplacianrecursivegaussianimagefilternew.h
)

SET( sd2d_SRCS
	seed_detection/seedsDetection_2D_Cpp.cxx
	seed_detection/seedsdetection.h
)

SET( mcgl3d_SRCS
	graphColLearn_3D/Multi_Color_Graph_Learning_3D.h
	graphColLearn_3D/gvc.h
	graphColLearn_3D/Multi_Color_Graph_Learning_3D.cpp
	graphColLearn_3D/sequential_coloring.cpp
)

SET( emgmm_SRCS
	EM_GMM/EM_Project_3_cpp_3D_comp.cxx
	EM_GMM/EM_Project_3_cpp_3D_comp.h
)

SET( clust_SRCS
	clustering_3D/local_max_clust_3D.h
	clustering_3D/Local_Max_Clust_3D_CPP.cpp
)

SET( border_SRCS
  cellLabelToBorders/label2border.h
  cellLabelToBorders/itkLabelBorderImageFilter.h
  cellLabelToBorders/label2border.cpp
  cellLabelToBorders/itkLabelBorderImageFilter.txx
)

SET( cbin_SRCS
  cell_binarization/cell_binarization.h
  cell_binarization/cell_binarization.cpp    
  cell_binarization/new_graph.h
  cell_binarization/block_b.h   
)

SET( cb3_SRCS
	alpha_expansion_3d/alpha_expansion_3d.h
	alpha_expansion_3d/alpha_expansion_3d.cpp
	alpha_expansion_3d/GCoptimization.h
	alpha_expansion_3d/GCoptimization.cpp	
	alpha_expansion_3d/LinkedBlockList.h
	alpha_expansion_3d/LinkedBlockList.cpp
	#alpha_expansion_3d/GraphCut3dConstr.cpp	
	#alpha_expansion_3d/GraphCutMex.cpp
	alpha_expansion_3d/graph.cpp
	alpha_expansion_3d/maxflow.cpp
	alpha_expansion_3d/block.h
	alpha_expansion_3d/graph.h
	alpha_expansion_3d/GraphCut.h
)

ADD_LIBRARY( Yousef_Nucleus_Seg ${NSRCS} ${sd3d_SRCS} ${sd2d_SRCS} ${mcgl3d_SRCS} ${emgmm_SRCS} ${clust_SRCS} ${border_SRCS} ${cbin_SRCS} ${cb3_SRCS})

TARGET_LINK_LIBRARIES ( Yousef_Nucleus_Seg ITKBasicFilters ITKCommon ITKIO)
