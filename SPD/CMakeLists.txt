cmake_minimum_required(VERSION 2.8)
PROJECT(SPD)

INCLUDE_DIRECTORIES(${VXL_MUL_INCLUDE_DIR}) 

SET( Libraries
    ${VXL_LIBRARIES}
    vtkCommon
    vtksys
    vtkFiltering
    vcl
    vnl
    mbl
    ${QT_LIBRARIES})

QT4_WRAP_CPP( GUI_MOC_SRCS spdmainwindow.h )

SET(SPD_SRCS
	SPDAnalysisModel.h
	SPDAnalysisModel.cpp
	transportSimplex.h)

SET(SPD_TEST_SRCS
	spdmainwindow.h
	spdmainwindow.cpp )

ADD_LIBRARY( SPD ${SPD_SRCS} )
TARGET_LINK_LIBRARIES(SPD ${Libraries} ftkGUI ftkCommon)
ADD_EXECUTABLE(SPDTEST main.cpp ${GUI_MOC_SRCS} ${SPD_TEST_SRCS} )
                
TARGET_LINK_LIBRARIES(SPDTEST SPD)


