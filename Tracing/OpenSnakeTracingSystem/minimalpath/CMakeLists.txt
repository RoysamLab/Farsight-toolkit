# The header files
SET(HEADERS
  itkSingleImageCostFunction.h
  itkArrivalFunctionToPathFilter.h
  itkSpeedFunctionToPathFilter.h
)

# The implementation files
SET(SOURCES
  itkSingleImageCostFunction.txx
  itkArrivalFunctionToPathFilter.txx
  itkSpeedFunctionToPathFilter.txx
  #examples.cxx
)

ADD_LIBRARY( MinimalPath
			${SOURCES}
			${HEADERS}
)
SET_TARGET_PROPERTIES(MinimalPath PROPERTIES LINKER_LANGUAGE C)