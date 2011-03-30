SUBDIRS( fuzzyClustering )
SUBDIRS( libsvm )
SUBDIRS( agf )

IF(BUILD_EMBREX)
  SUBDIRS( embrex )
ENDIF(BUILD_EMBREX)