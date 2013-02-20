PROJECT(Haralick)


ADD_LIBRARY( Haralick  Haralick.cpp Haralick.h )
ADD_EXECUTABLE(Haralick_test Haralick_test.cpp )

TARGET_LINK_LIBRARIES(Haralick_test Haralick             ${ITK_LIBRARIES} )




