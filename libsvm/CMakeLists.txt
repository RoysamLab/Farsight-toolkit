ADD_LIBRARY( libSVM svm.cpp svm.h )

ADD_EXECUTABLE(svm-predict svm-predict.c)
TARGET_LINK_LIBRARIES(svm-predict libSVM)
ADD_EXECUTABLE(svm-scale svm-scale.c)
TARGET_LINK_LIBRARIES(svm-scale libSVM)
ADD_EXECUTABLE(svm-train svm-train.c)
TARGET_LINK_LIBRARIES(svm-train libSVM)


INSTALL(TARGETS svm-predict svm-scale svm-train DESTINATION bin)