SET( GUI_SRCS
    	PlotWindow.cpp
    	ScatterView.cpp
    	TableWindow.cpp
	FTKRenderWindow.cpp
	HistoWindow.cpp
    	ImageBrowser5D.cpp
    	PatternAnalysisWizard.cpp
    	GenericALDialog.cpp
    	Active_Learning_Dialog.cpp
    	ObjectSelection.cpp
    	LabelImageViewQT.cpp
    	TrainingDialog.cpp
	ftkPreprocessDialog.cpp
	PreprocessDialog.cpp
	StatisticsToolbar.cpp
	GraphWindowForNewSelection.cpp
	GraphWindow.cpp
	SelectiveClustering.cpp
	QvtkTableView.cpp
	SelectionUtilities.cpp
)
SET( NONQT_GUI_SRCS
	ProcessObjectProgressUpdater.cpp
)

SET( GUI_HDRS
    	PlotWindow.h
    	ScatterView.h
    	TableWindow.h
	FTKRenderWindow.h
	HistoWindow.h
    	ImageBrowser5D.h
    	PatternAnalysisWizard.h
    	GenericALDialog.h
    	Active_Learning_Dialog.h
    	ObjectSelection.h
    	LabelImageViewQT.h
    	TrainingDialog.h
	ftkPreprocessDialog.h
	PreprocessDialog.h
	StatisticsToolbar.h
	GraphWindowForNewSelection.h
	GraphWindow.h
	SelectiveClustering.h
	QvtkTableView.h
	SelectionUtilities.h
)
SET( NONQT_GUI_HDRS
	ProcessObjectProgressUpdater.h
)

SET( GUI_LIBS
  ${QT_LIBRARIES} QVTK ftkCommon ftkPreprocess ftkImage TinyXML libSVM
  vtkVolumeRendering vtkHybrid vtkCharts vnl)

IF(QtTestingFound)
  CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/GUITesterHelper.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/GUITesterHelper.h
  @ONLY)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
  
  # normally this directory is created automatically by CTest, but we do it
  # manually here for the sake of nightly superbuilds
  IF(NOT EXISTS ${CMAKE_BINARY_DIR}/Testing/Temporary)
    FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary)
  ENDIF()
  
  LIST( APPEND GUI_SRCS GUITester.cpp )
  LIST( APPEND GUI_HDRS GUITester.h )
  LIST( APPEND GUI_LIBS QtTesting )
ENDIF()

IF(USE_KPLS)
  LIST( APPEND GUI_LIBS ${LIB_KPLS} )
ENDIF(USE_KPLS)

SET( GUI_UIS
)

SET( GUI_RCS
)

# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
QT4_ADD_RESOURCES( GUI_RC_SRCS ${GUI_RCS} )

# this will run uic on .ui files:
QT4_WRAP_UI( GUI_UI_HDRS ${GUI_UIS} )

# and finally this will run moc:
QT4_WRAP_CPP( GUI_MOC_SRCS ${GUI_HDRS} )

# This must be built static because libKPLS is only available as a static
# library.
ADD_LIBRARY( ftkGUI STATIC
  ${GUI_SRCS}
  ${GUI_MOC_SRCS}
  ${GUI_HDRS}
  ${GUI_UI_HDRS}
  ${GUI_RC_SRCS}
  ${NONQT_GUI_SRCS}
  ${NONQT_GUI_HDRS}
  )
# For linking into shared libraries.
IF( CMAKE_COMPILER_IS_GNUCC )
  SET_TARGET_PROPERTIES( ftkGUI PROPERTIES COMPILE_FLAGS -fPIC )
ENDIF()

TARGET_LINK_LIBRARIES( ftkGUI ${GUI_LIBS} )
