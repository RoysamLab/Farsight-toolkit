project(shtensor)

add_executable(scalar_voting scalar_voting.cpp)
target_link_libraries(scalar_voting ITKIO ITKCommon)
install(TARGETS scalar_voting RUNTIME DESTINATION bin)
