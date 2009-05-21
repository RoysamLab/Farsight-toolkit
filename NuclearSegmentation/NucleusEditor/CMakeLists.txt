SET( NucluesEditor_SRCS
	NucleusEditor.cpp
    	SegmentationView.cpp
    	SegmentationWindow.cpp
    	SegmentationModel.cpp
   	NuclearSegmentationWizard.cpp
)

SET( NucluesEditor_HDRS
	NucleusEditor.h
    	SegmentationView.h
    	SegmentationWindow.h
    	SegmentationModel.h
    	NuclearSegmentationWizard.h
)

SET( NucluesEditor_UIS
)

SET( NucluesEditor_RCS
)

# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
QT4_ADD_RESOURCES( NucluesEditor_RC_SRCS ${NucluesEditor_RCS} )

# this will run uic on .ui files:
QT4_WRAP_UI( NucluesEditor_UI_HDRS ${NucluesEditor_UIS} )

# and finally this will run moc:
QT4_WRAP_CPP( NucluesEditor_MOC_SRCS ${NucluesEditor_HDRS} )

#ADD_LIBRARY( NucluesEditor STATIC 
#             ${NucluesEditor_SRCS}
#             ${NucluesEditor_HDRS} 
#             ${NucluesEditor_MOC_SRCS} 
#             ${NucluesEditor_UI_HDRS}
#             ${NucluesEditor_RC_SRCS}
#)

ADD_EXECUTABLE(NucleusEditor main.cpp
             ${NucluesEditor_SRCS}
             ${NucluesEditor_HDRS} 
             ${NucluesEditor_MOC_SRCS} 
             ${NucluesEditor_UI_HDRS}
             ${NucluesEditor_RC_SRCS}
)

TARGET_LINK_LIBRARIES( NucleusEditor ${QT_LIBRARIES} QVTK ftkImage TinyXML libSVM vtkVolumeRendering NuclearSegmentation ftkGUI)

#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
SET_TARGET_PROPERTIES( NuclearSegmentation PROPERTIES WIN32_EXECUTABLE 1)
#OR MAKE A CONSOLE APPLICATION
#SET_TARGET_PROPERTIES( NuclearSegmentation PROPERTIES WIN32_EXECUTABLE 0)

INSTALL( TARGETS NucleusEditor DESTINATION bin )


