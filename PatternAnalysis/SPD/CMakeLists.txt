SET( Libraries
    ${VXL_LIBRARIES}
    ${VTK_LIBRARIES}
    ${ITK_LIBRARIES}
    vcl
    vnl
    mbl
    ${QT_LIBRARIES})

QT4_WRAP_CPP( GUI_MOC_SRCS SPDMSTModuleMatch.h SPDkNNGModuleMatch.h HeatmapWindow.h HeatmapWindowForNewSelection.h ProgressionHeatmapWindow.h)

SET(SPD_SRCS
	SPDAnalysisModel.h
	SPDAnalysisModel.cpp
	transportSimplex.h
	HeatmapWindow.h
	HeatmapWindow.cpp
	HeatmapWindowForNewSelection.h
	HeatmapWindowForNewSelection.cpp
	ProgressionHeatmapWindow.h
	ProgressionHeatmapWindow.cpp
	SPDkNNGModuleMatch.h
	SPDkNNGModuleMatch.cpp
	SPDMSTModuleMatch.h
	SPDMSTModuleMatch.cpp)

IF ( NOT BUILD_CLUSCLUS )
	MESSAGE( SEND_ERROR "Must turn on BUILD_CLUSCLUS" )
ENDIF ( NOT BUILD_CLUSCLUS )
	
ADD_LIBRARY( SPD ${SPD_SRCS} ${GUI_MOC_SRCS})
TARGET_LINK_LIBRARIES(SPD ${Libraries} ftkGUI ftkCommon ClusClus_lib)
ADD_EXECUTABLE(PSC SPDAnalysisModel.cpp SPDAnalysisModel.h transportSimplex.h PSC.cpp)
TARGET_LINK_LIBRARIES(PSC vnl mbl vul ${VTK_LIBRARIES} ftkCommon ClusClus_lib)