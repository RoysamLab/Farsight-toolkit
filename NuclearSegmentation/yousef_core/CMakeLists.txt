SET( NSRCS
	yousef_seg.h
	yousef_seg.cpp
  gvc/gvc.h
  gvc/gvc.cpp
)

################################################################################
#      Here we add all sub-directories that are involved with this program       #
################################################################################

SET( sd3d_SRCS
	seed_detection_3D/seedsDetection_3D.cxx
	seed_detection_3D/seedsdetection_3D.h
)

SET( sd2d_SRCS
	seed_detection/seedsDetection_2D_Cpp.cxx
	seed_detection/seedsdetection.h
)

SET( mcgl3d_SRCS
	graphColLearn_3D/Multi_Color_Graph_Learning_3D.h
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
SET( clust2d_SRCS
	clustering_2D/Local_Max_Clust_2D.h
	clustering_2D/Local_Max_Clust_2D.cpp
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
  cell_binarization/adaptive_binarization.h
  cell_binarization/adaptive_binarization.cpp
  cell_binarization/new_graph.h
  cell_binarization/block_b.h   
)

SET( cb3_SRCS
	alpha_expansion/alpha_expansion.h
	alpha_expansion/alpha_expansion.cpp
	alpha_expansion/GCoptimization.h
	alpha_expansion/GCoptimization.cpp	
	alpha_expansion/LinkedBlockList.h
	alpha_expansion/LinkedBlockList.cpp
	alpha_expansion/graph.cpp
	alpha_expansion/maxflow.cpp
	alpha_expansion/block.h
	alpha_expansion/graph.h
	alpha_expansion/GraphCut.h
	alpha_expansion/Multi_Color_Graph_Learning_2D.h
	alpha_expansion/Multi_Color_Graph_Learning_2d.cpp
)

SET( prepr_SRCS
	preprocessing/GrAnisDiff.h
	preprocessing/GrAnisDiff.cpp
)

INCLUDE_DIRECTORIES(gvc)


#######################################
# ADD CUDA FILES HERE
#######################################


if (BUILD_CUDA)
    SET ( sd3d_SRCS ${sd3d_SRCS} seed_detection_3D/LocalMaxima3D.cu  )
    SET( clust_SRCS ${clust_SRCS} clustering_3D/Local_Max_Clust_3D.cu )

    CUDA_ADD_LIBRARY( Yousef_Nucleus_Seg ${NSRCS} ${sd3d_SRCS} ${sd2d_SRCS} ${mcgl3d_SRCS} ${emgmm_SRCS} ${clust_SRCS} ${clust2d_SRCS} ${border_SRCS} ${cbin_SRCS} ${cb3_SRCS} ${prepr_SRCS})
else (BUILD_CUDA)
    ADD_LIBRARY( Yousef_Nucleus_Seg ${NSRCS} ${sd3d_SRCS} ${sd2d_SRCS} ${mcgl3d_SRCS} ${emgmm_SRCS} ${clust_SRCS} ${clust2d_SRCS} ${border_SRCS} ${cbin_SRCS} ${cb3_SRCS} ${prepr_SRCS})
endif (BUILD_CUDA)

TARGET_LINK_LIBRARIES ( Yousef_Nucleus_Seg    ftkCommon     ${ITK_LIBRARIES} )
