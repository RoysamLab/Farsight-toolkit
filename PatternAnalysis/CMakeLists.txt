SUBDIRS( fuzzyClustering )
SUBDIRS( libsvm )

IF(BUILD_EMBREX)
  SUBDIRS( embrex )
ENDIF(BUILD_EMBREX)