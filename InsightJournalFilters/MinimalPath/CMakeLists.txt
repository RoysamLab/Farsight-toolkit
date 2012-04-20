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
  main.cxx
  #examples.cxx
)

# Add this as include directory
INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}
  ${SOURCE_PATH}
  ${VXL_INCLUDE_DIRS}
)

# Main library
ADD_EXECUTABLE(MinimalPathMain ${HEADERS} ${SOURCES})
TARGET_LINK_LIBRARIES(MinimalPathMain ${ITK_LIBRARIES})
