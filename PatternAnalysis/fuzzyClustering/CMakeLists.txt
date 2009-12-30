ADD_EXECUTABLE(run_fuzzy_clustering run_fuzzy_clustering.cpp fuzzy_clustering.cpp fuzzy_clustering.h)

INSTALL(TARGETS run_fuzzy_clustering DESTINATION ${INSTALL_BIN_DIR} )
