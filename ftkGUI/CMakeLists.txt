SET( GUI_SRCS
    	PlotWindow.cpp
    	ScatterView.cpp
    	TableWindow.cpp	
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
	GraphWindow.cpp
	
)

SET( GUI_HDRS
    	PlotWindow.h
    	ScatterView.h
    	TableWindow.h
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
	GraphWindow.h
	ColorMap.h
)


IF(BUILD_EMBREX)
SET( LIB_KPLS libKPLS )
ENDIF(BUILD_EMBREX)

SET( GUI_UIS
)

SET( GUI_RCS
)

#SET(QT_USE_OPENGL TRUE)

# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
QT4_ADD_RESOURCES( GUI_RC_SRCS ${GUI_RCS} )

# this will run uic on .ui files:
QT4_WRAP_UI( GUI_UI_HDRS ${GUI_UIS} )

# and finally this will run moc:
QT4_WRAP_CPP( GUI_MOC_SRCS ${GUI_HDRS} )

ADD_LIBRARY( ftkGUI
             ${GUI_SRCS}
             ${GUI_HDRS}
             ${GUI_MOC_SRCS}
             ${GUI_UI_HDRS}
             ${GUI_RC_SRCS}
)

TARGET_LINK_LIBRARIES( ftkGUI
  ${QT_LIBRARIES} QVTK ftkCommon ftkPreprocess ftkImage TinyXML libSVM vtkVolumeRendering vtkHybrid vnl
) 

IF(USE_KPLS)
  TARGET_LINK_LIBRARIES( ftkGUI ${LIB_KPLS} )
ENDIF(USE_KPLS)
