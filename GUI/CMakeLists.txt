INCLUDE_DIRECTORIES(
  ${PROJECT_BINARY_DIR}/GUI
#  ${PROJECT_BINARY_DIR}/Qt/Python
#  ${PROJECT_SOURCE_DIR}/Qt/Python
)

SET( GUI_SRCS
    ControlBar.cpp
    PlotWindow.cpp
    ScatterView.cpp
    TableWindow.cpp
    SegmentationView.cpp
    SegmentationWindow.cpp
    ModuleWidget.cpp
    SegmentationModel.cpp
    
    #ObjectSelection.cpp
    ImageBrowser5D.cpp
)

SET( GUI_HDRS
    ControlBar.h
    PlotWindow.h
    ScatterView.h
    TableWindow.h
    SegmentationView.h
    SegmentationWindow.h
    ModuleWidget.h
    SegmentationModel.h
    
    #ObjectSelection.h
    ImageBrowser5D.h
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

ADD_LIBRARY( GUI STATIC 
             ${GUI_SRCS}
             ${GUI_HDRS} 
             ${GUI_MOC_SRCS} 
             ${GUI_UI_HDRS}
             ${GUI_RC_SRCS}
)

TARGET_LINK_LIBRARIES( GUI ${QT_LIBRARIES} QVTK ftkImage TinyXML libSVM vtkVolumeRendering)
#TARGET_LINK_LIBRARIES( GUI ${QT_LIBRARIES} FTKImage TinyXML QtPython)

#########################################
# ADD NEW MODULES HERE
# FOLLOW BUILD_YOUSEF EXAMPLE
#########################################

IF(BUILD_NUCLEI)
	TARGET_LINK_LIBRARIES( GUI NuclearSegmentation )
ENDIF(BUILD_NUCLEI)	

############################################3


# last thing we have to do is to tell CMake what libraries our executable needs,
ADD_EXECUTABLE( Farsight main.cpp )

#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
SET_TARGET_PROPERTIES(Farsight PROPERTIES WIN32_EXECUTABLE 1)
#OR MAKE A CONSOLE APPLICATION
#SET_TARGET_PROPERTIES(Farsight PROPERTIES WIN32_EXECUTABLE 0)

TARGET_LINK_LIBRARIES( Farsight 
                       GUI
)

INSTALL( TARGETS Farsight DESTINATION bin )
#INSTALL( TARGETS Farsight GUI DESTINATION . )

