cmake_minimum_required(VERSION 2.8)
PROJECT(SPD)

INCLUDE_DIRECTORIES(${VXL_MUL_INCLUDE_DIR}) 

SET( Libraries
    ${VXL_LIBRARIES}
    ${VTK_LIBRARIES}
    vcl
    vnl
    mbl
    ${QT_LIBRARIES})

QT4_WRAP_CPP( GUI_MOC_SRCS spdtestwindow.h spdtestwindowForNewSelection.h HeatmapWindow.h HeatmapWindowForNewSelection.h ProgressionHeatmapWindow.h)

SET(SPD_SRCS
	SPDAnalysisModel.h
	SPDAnalysisModel.cpp
	transportSimplex.h
	LCSAnalysis.h
	LCSAnalysis.cpp
	spdtestwindow.h
	spdtestwindow.cpp
	spdtestwindowForNewSelection.h
	spdtestwindowForNewSelection.cpp
	HeatmapWindow.h
	HeatmapWindow.cpp
	HeatmapWindowForNewSelection.h
	HeatmapWindowForNewSelection.cpp
	ProgressionHeatmapWindow.h
	ProgressionHeatmapWindow.cpp)

IF ( NOT BUILD_CLUSCLUS )
	MESSAGE( SEND_ERROR "Must turn on BUILD_CLUSCLUS" )
ENDIF ( NOT BUILD_CLUSCLUS )
	
ADD_LIBRARY( SPD ${SPD_SRCS} ${GUI_MOC_SRCS})
TARGET_LINK_LIBRARIES(SPD ${Libraries} ftkGUI ftkCommon ClusClus_lib)
ADD_EXECUTABLE(PSC SPDAnalysisModel.cpp SPDAnalysisModel.h transportSimplex.h PSC.cpp)
TARGET_LINK_LIBRARIES(PSC vnl mbl ${VTK_LIBRARIES} ftkCommon ClusClus_lib)