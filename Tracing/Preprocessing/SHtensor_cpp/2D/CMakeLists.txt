project(shtensor)

include(FindFFTW3.cmake)
include_directories(${FFTW3_INCLUDE_DIR})
link_directories(${FFTW3_LINK_DIRECTORIES})

add_executable(scalar_voting_2d scalar_voting_2d.cpp)
add_executable(tensor_voting_2d tensor_voting_2d.cpp)
target_link_libraries(scalar_voting_2d ITKIO ITKCommon vnl vnl_algo ${FFTW3_LIBRARIES})
target_link_libraries(tensor_voting_2d ITKIO ITKCommon vnl vnl_algo ${FFTW3_LIBRARIES})
install(TARGETS scalar_voting_2d RUNTIME DESTINATION bin)
install(TARGETS tensor_voting_2d RUNTIME DESTINATION bin)
