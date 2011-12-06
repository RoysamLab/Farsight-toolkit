SET( SampleEditor_SRCS
	SampleEditor.cpp
)

SET( SampleEditor_HDRS
	SampleEditor.h
)

SET( SampleEditor_UIS
)

SET( SampleEditor_RCS
)

# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
QT4_ADD_RESOURCES( SampleEditor_RC_SRCS ${SampleEditor_RCS} )

# this will run uic on .ui files:
QT4_WRAP_UI( SampleEditor_UI_HDRS ${SampleEditor_UIS} )

# and finally this will run moc:
QT4_WRAP_CPP( SampleEditor_MOC_SRCS ${SampleEditor_HDRS} )


ADD_EXECUTABLE(SampleEditor main.cpp
             ${SampleEditor_SRCS}
             ${SampleEditor_HDRS} 
             ${SampleEditor_MOC_SRCS} 
             ${SampleEditor_UI_HDRS}
             ${SampleEditor_RC_SRCS}
)

TARGET_LINK_LIBRARIES( SampleEditor ${QT_LIBRARIES} ftkCommon ftkGraphs SPD ClusClus_lib)

#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
#SET_TARGET_PROPERTIES( SampleEditor PROPERTIES WIN32_EXECUTABLE 1)
#OR MAKE A CONSOLE APPLICATION
SET_TARGET_PROPERTIES( SampleEditor PROPERTIES WIN32_EXECUTABLE 0)

INSTALL( TARGETS SampleEditor DESTINATION ${INSTALL_BIN_DIR} )


