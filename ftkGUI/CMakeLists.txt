SET( GUI_SRCS
    	#ControlBar.cpp
    	PlotWindow.cpp
    	ScatterView.cpp
    	TableWindow.cpp		
   	#TableView.cpp
   	#Histogram.cpp
   		HistoWindow.cpp
    	ImageBrowser5D.cpp
    	LibSVMWidget.cpp
    	PatternAnalysisWizard.cpp
)

SET( GUI_HDRS
    	#ControlBar.h
    	PlotWindow.h
    	ScatterView.h
    	TableWindow.h
   		#TableView.h
    	#Histogram.h
    	HistoWindow.h
    	ImageBrowser5D.h
    	LibSVMWidget.h
    	PatternAnalysisWizard.h
)

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

ADD_LIBRARY( ftkGUI STATIC 
             ${GUI_SRCS}
             ${GUI_HDRS} 
             ${GUI_MOC_SRCS} 
             ${GUI_UI_HDRS}
             ${GUI_RC_SRCS}
)

TARGET_LINK_LIBRARIES( ftkGUI ${QT_LIBRARIES} QVTK ftkImage TinyXML libSVM libKPLS vtkVolumeRendering)

