SET( NucleusEditor_SRCS
  NucleusEditor.cpp
  ProjectFilenamesDialog.cpp
  ExclusionDialog.cpp
  PreferencesDialog.cpp
  ftkProjectDefinition.cpp
  ftkProjectProcessor.cpp
  ftkProjectFiles.cpp
  ActiveLearningDialog.cpp
  GalleryDialog.cpp	
)

SET( NucleusEditor_HDRS
  NucleusEditor.h
  ProjectFilenamesDialog.h
  ExclusionDialog.h
  PreferencesDialog.h
  ftkProjectDefinition.h
  ftkProjectProcessor.h
  ftkProjectFiles.h
  ActiveLearningDialog.h
  GalleryDialog.h
)

SET( NucleusEditor_MOC_REQD_HDRS
  NucleusEditor.h 
  ProjectFilenamesDialog.h 
  ExclusionDialog.h
  PreferencesDialog.h
  ActiveLearningDialog.h
  GalleryDialog.h  
)

SET( NucleusEditor_UIS
)

SET( NucleusEditor_RCS
  icons.qrc
)

IF(BUILD_CLUSCLUS)
SET (CLUSCLUS_LIBRARIES
    ClusClus_lib
    )
ENDIF (BUILD_CLUSCLUS)

# This part makes sure to include the necessary modification 
# if tracking is to be used with nucleus editor:
IF (USE_TRACKING)
SET( NucleusEditor_SRCS
	NucleusEditor.cpp
	ProjectFilenamesDialog.cpp
	ExclusionDialog.cpp
	PreferencesDialog.cpp
	ftkProjectDefinition.cpp
	ftkProjectProcessor.cpp
	ftkProjectFiles.cpp
	ActiveLearningDialog.cpp
	GalleryDialog.cpp	
	KymoGraphView.cpp
	Trace.cpp
	TrackingDialog.cpp
	Image3DView.cpp
	)		
SET( NucleusEditor_HDRS
    NucleusEditor.h
    ProjectFilenamesDialog.h
    ExclusionDialog.h
    PreferencesDialog.h
    ftkProjectDefinition.h
    ftkProjectProcessor.h
    ftkProjectFiles.h
    ActiveLearningDialog.h
    GalleryDialog.h
	KymoGraphView.h
	Trace.h
	TrackingDialog.h
	Image3DView.h
	)
SET( NucleusEditor_MOC_REQD_HDRS
    NucleusEditor.h 
    ProjectFilenamesDialog.h 
    ExclusionDialog.h
    PreferencesDialog.h
    ActiveLearningDialog.h
    GalleryDialog.h  
	KymoGraphView.h
	TrackingDialog.h
	Image3DView.h
	)
	FIND_PACKAGE(CPLEX)
	FIND_PACKAGE(CONCERT)
	IF(CPLEX_FOUND)
		INCLUDE_DIRECTORIES(${CPLEX_INCLUDE_DIRS})
	ENDIF(CPLEX_FOUND)
	IF(CONCERT_FOUND)
		INCLUDE_DIRECTORIES(${CONCERT_INCLUDE_DIRS})
	ENDIF(CONCERT_FOUND)
	ADD_DEFINITIONS(-DIL_STD)
	
ENDIF (USE_TRACKING)

# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
QT4_ADD_RESOURCES( NucleusEditor_RC_SRCS ${NucleusEditor_RCS} )

# this will run uic on .ui files:
QT4_WRAP_UI( NucleusEditor_UI_HDRS ${NucleusEditor_UIS} )

# and finally this will run moc:
QT4_WRAP_CPP( NucleusEditor_MOC_SRCS ${NucleusEditor_MOC_REQD_HDRS} )

ADD_EXECUTABLE( NucleusEditor ${EXE_TYPE}
                main.cpp
                ${NucleusEditor_SRCS}
                ${NucleusEditor_HDRS} 
                ${NucleusEditor_MOC_SRCS} 
                ${NucleusEditor_UI_HDRS}
                ${NucleusEditor_RC_SRCS} )


SET (NUCLEUSEDITOR_LIBRARIES
	${ITK_USE_LIBRARIES}
	${QT_LIBRARIES}
	ftkGraphs
	ftkGUI 
	ftkCommon 
	ftkFeatures 
	ftkSpectralUnmixing
	libSVM 
	vtkVolumeRendering 
	QVTK 
	NuclearSegmentation 
	CytoplasmSegmentation 
	PixelLevelAnalysis 
	NESQLITE 
	sqlite3 
	MCLR
	${CLUSCLUS_LIBRARIES}
	)
IF (USE_TRACKING)
	TARGET_LINK_LIBRARIES( NucleusEditor ${NUCLEUSEDITOR_LIBRARIES} CellTrackerLib)
ELSE (USE_TRACKING)
	TARGET_LINK_LIBRARIES( NucleusEditor ${NUCLEUSEDITOR_LIBRARIES} ftkImage)
ENDIF (USE_TRACKING)

ADD_LIBRARY( Project_Processor ftkProjectProcessor.cpp ftkProjectProcessor.h ftkProjectDefinition.cpp ftkProjectDefinition.h)
TARGET_LINK_LIBRARIES( Project_Processor ${QT_LIBRARIES}  ftkGraphs ftkGUI ftkImage ftkCommon ftkFeatures NuclearSegmentation CytoplasmSegmentation PixelLevelAnalysis NESQLITE sqlite3 MCLR)

ADD_EXECUTABLE( projproc processor_main.cpp ftkProjectProcessor.cpp ftkProjectDefinition.cpp)
TARGET_LINK_LIBRARIES( projproc ${QT_LIBRARIES}  ftkGraphs ftkGUI ftkImage ftkCommon ftkFeatures NuclearSegmentation CytoplasmSegmentation PixelLevelAnalysis NESQLITE sqlite3 MCLR)
install(TARGETS projproc DESTINATION bin )

IF (BUILD_MODEL_SEG)
    IF (MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D MODEL_SEG")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D MODEL_SEG")
	ENDIF(MSVC)
	TARGET_LINK_LIBRARIES( NucleusEditor ModelSeg )
	TARGET_LINK_LIBRARIES( projproc ModelSeg )
ENDIF(BUILD_MODEL_SEG)

#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
#SET_TARGET_PROPERTIES( NucleusEditor PROPERTIES WIN32_EXECUTABLE 1)
#OR MAKE A CONSOLE APPLICATION
SET_TARGET_PROPERTIES( NucleusEditor PROPERTIES WIN32_EXECUTABLE 0)

if(APPLE)
   install(TARGETS NucleusEditor DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS NucleusEditor DESTINATION bin )
endif(APPLE)

