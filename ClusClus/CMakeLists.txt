cmake_minimum_required(VERSION 2.8)
PROJECT(ClusClus)

SET( Libraries
    vtkCommon
    vtksys
    vtkFiltering
    vtkVolumeRendering
    vtkViews
    ${QT_LIBRARIES})

SET( ClusClus_SRCS  
	Dendrogram.cpp
	clusclus.cpp
	clusgap.cpp
	clusclusMainwindow.cpp
	HeatmapWindow.cpp
	ProgressionHeatmapWindow.cpp)

SET( ClusClus_HDRS 
	Dendrogram.h
	clusclus.h
	clusgap.h
	clusclusMainwindow.h
	HeatmapWindow.h
	ProgressionHeatmapWindow.h)
	
LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS})
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS})
QT4_WRAP_CPP( ClusClus_MOC_REQD_HDRS HeatmapWindow.h Dendrogram.h clusclusMainwindow.h ProgressionHeatmapWindow.h)  

ADD_LIBRARY( ClusClus_lib ${ClusClus_SRCS} ${ClusClus_HDRS} ${ClusClus_MOC_REQD_HDRS} )

TARGET_LINK_LIBRARIES( ClusClus_lib  
    ${Libraries}
    ftkGUI
    ftkCommon)           

ADD_EXECUTABLE( ClusClus main.cpp )
  
TARGET_LINK_LIBRARIES( ClusClus ClusClus_lib)

