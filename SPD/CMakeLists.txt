cmake_minimum_required(VERSION 2.8)
PROJECT(SPD_CMAKE)

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

QT4_WRAP_CPP( GUI_MOC_SRCS SPDMainWindow.h)

SET( SPD_SRCS
	spdmainwindow.cpp
	SPDAnalysisModel.cpp
)

SET( SPD_HDRS
	spdmainwindow.h
	SPDAnalysisModel.h
)

ADD_EXECUTABLE(SPD main.cpp 
                ${SPD_SRCS}
				${SPD_HDRS}
                ${GUI_MOC_SRCS})
                
TARGET_LINK_LIBRARIES(SPD ${Libraries} ftkGUI)


