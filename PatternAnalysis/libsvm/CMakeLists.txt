IF(UNIX)
  SET(MY_MATH_LIB -lm)
ELSE(UNIX)
  SET(MY_MATH_LIB)
ENDIF(UNIX)

ADD_LIBRARY( libSVM svm.cpp svm.h)

#ADD_EXECUTABLE(svm-predict svm-predict.cpp)
#TARGET_LINK_LIBRARIES(svm-predict libSVM ${MY_MATH_LIB})
#ADD_EXECUTABLE(svm-scale svm-scale.cpp)
#TARGET_LINK_LIBRARIES(svm-scale libSVM ${MY_MATH_LIB})
#ADD_EXECUTABLE(svm-train svm-train.cpp)
#TARGET_LINK_LIBRARIES(svm-train libSVM ${MY_MATH_LIB})

#INSTALL(TARGETS svm-predict svm-scale svm-train DESTINATION ${INSTALL_BIN_DIR})
