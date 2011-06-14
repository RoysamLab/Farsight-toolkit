SUBDIRS( fuzzyClustering )
SUBDIRS( libsvm )
SUBDIRS( activeLearning )
SUBDIRS( agf )

IF(BUILD_EMBREX)
  SUBDIRS( embrex )
ENDIF(BUILD_EMBREX)