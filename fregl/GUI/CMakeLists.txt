INCLUDE_DIRECTORIES( ${Farsight_BINARY_DIR}/fregl/GUI )

SET( MY_GUI_SRCS
    gnrcImageDialog.cpp
    MainWindow.cpp
    deluxeScene.cpp
)

SET( MY_GUI_HDRS
    gnrcimagedialog.h
    mainwindow.h
    deluxescene.h
)


SET( MY_UIS
   ImageViewArea.ui
)


SET(QT_USE_OPENGL TRUE)

# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
#QT4_ADD_RESOURCES( MY_RC_SRCS ${MY_RCS} )

# this will run uic on .ui files:
QT4_WRAP_UI( MY_UI_HDRS ${MY_UIS} )

# and finally this will run moc:
QT4_WRAP_CPP( MY_MOC_SRCS ${MY_GUI_HDRS} )

ADD_LIBRARY( fregl_gui
             ${MY_GUI_SRCS}
	     ${MY_GUI_HDRS}	 
             ${MY_MOC_SRCS} 
             ${MY_UI_HDRS}
             #${MY_RC_SRCS}
)

TARGET_LINK_LIBRARIES( fregl_gui ${QT_LIBRARIES} fregl)

#This line insures that the window is opened as a standalone function. If we 
#need a console then we add -mconsole after -mwindows
#set(LINK_FLAGS "${LINK_FLAGS} -mwindows -mconsole")

#LINK_LIBRARIES(
#	vtksys
#    	vtkCommon
#    	vtkIO
#    	vtkFiltering
#)

# last thing we have to do is to tell CMake what libraries our executable needs,
# luckily FIND_PACKAGE prepared QT_LIBRARIES variable for us:

ADD_EXECUTABLE( MontageNavigator main.cpp)
TARGET_LINK_LIBRARIES( MontageNavigator 
                       #${LINK_FLAGS} 
                       fregl_gui
)

INSTALL( TARGETS MontageNavigator DESTINATION ${INSTALL_BIN_DIR} )
