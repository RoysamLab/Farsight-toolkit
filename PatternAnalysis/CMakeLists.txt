SUBDIRS( fuzzyClustering )
SUBDIRS( libsvm )
SUBDIRS( activeLearning )

IF(BUILD_EMBREX)
  SUBDIRS( embrex )
ENDIF(BUILD_EMBREX)