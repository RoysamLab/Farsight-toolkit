project(shtensor)

include(FindFFTW3.cmake)
include_directories(${FFTW3_INCLUDE_DIR})
link_directories(${FFTW3_LINK_DIRECTORIES})
IF(WIN32)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
ELSE(WIN32)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF(WIN32)

add_executable(scalar_voting_2d scalar_voting_2d.cpp)
add_executable(tensor_voting_2d tensor_voting_2d.cpp)

target_link_libraries(scalar_voting_2d ${ITK_LIBRARIES} vnl vnl_algo ${FFTW3_LIBRARIES})
target_link_libraries(tensor_voting_2d ${ITK_LIBRARIES} vnl vnl_algo ${FFTW3_LIBRARIES})
install(TARGETS scalar_voting_2d RUNTIME DESTINATION ${INSTALL_BIN_DIR})
install(TARGETS tensor_voting_2d RUNTIME DESTINATION ${INSTALL_BIN_DIR})
