add_executable(MedianImageFilter MedianImageFilter.cxx)
target_link_libraries(MedianImageFilter ${ITK_LIBRARIES})

add_executable(RobustAutomaticThreshold RobustAutomaticThreshold.cxx)
target_link_libraries(RobustAutomaticThreshold ${ITK_LIBRARIES})

add_executable(HoleFilling HoleFilling.cxx)
target_link_libraries(HoleFilling ${ITK_LIBRARIES})

add_executable(Open Open.cxx)
target_link_libraries(Open ${ITK_LIBRARIES})

add_executable(IsolateCentroids IsolateCentroids.cxx)
target_link_libraries(IsolateCentroids ${ITK_LIBRARIES})
