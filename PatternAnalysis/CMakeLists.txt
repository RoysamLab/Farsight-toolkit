add_subdirectory( fuzzyClustering )
add_subdirectory( libsvm )
add_subdirectory( activeLearning )
add_subdirectory( agf )
add_subdirectory( DiffusionMap )

option( BUILD_GT_CLUSTERING "Build GT clustering" OFF )
if(BUILD_GT_CLUSTERING)
  set(CPACK_PACKAGE_EXECUTABLES "ftkGTClustering" "GT Clustering" ${CPACK_PACKAGE_EXECUTABLES})
  add_subdirectory( ftkGTClustering )
  add_definitions( -DUSE_GT_CLUSTERING )
  include_directories( ${PROJECT_SOURCE_DIR}/PatternAnalysis/ftkGTClustering )
endif(BUILD_GT_CLUSTERING)

