PROJECT(QtImageCombiner)

#Build the QtImageCombiner class files

ADD_LIBRARY(QtImageCombiner QtImageCombiner.cpp  QtImageCombiner.h )

ADD_EXECUTABLE(QtImageCombiner_main  QtImageCombiner_main.cpp  )

TARGET_LINK_LIBRARIES(QtImageCombiner_main  QtImageCombiner   ${QT_LIBRARIES}) 



