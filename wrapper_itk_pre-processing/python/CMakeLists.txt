PROJECT(wrapper_itk_pre-processing)

OPTION(USE_WRAPITK "Use WrapITK for Python wrappers; if OFF, the wrappers will be built using CSwig only." ON)    
INSTALL(DIRECTORY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/wrapper_itk_pre-processing/python FILES_MATCHING PATTERN "*.py" PATTERN "ExecuteAlgorithm*.py" EXCLUDE)
INSTALL(DIRECTORY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/wrapper_itk_pre-processing/python FILES_MATCHING PATTERN "*.xml")

IF(USE_WRAPITK)
    INSTALL(FILES ExecuteAlgorithmModule.py ExecuteAlgorithmWrapITK.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/wrapper_itk_pre-processing/python)
ELSE(USE_WRAPITK)
    INSTALL(FILES ExecuteAlgorithmModule.py ExecuteAlgorithmCSwig.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/wrapper_itk_pre-processing/python)
ENDIF(USE_WRAPITK)    
