SET( HistoGUI_SRCS
  HistoGUI.cpp
  ftkHistopathology.cpp
)

SET( HistoGUI_HDRS
  HistoGUI.h
  ftkHistopathology.h
)

SET( HistoGUI_UIS
)

SET( HistoGUI_RCS
)

# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
QT4_ADD_RESOURCES( HistoGUI_RC_SRCS ${HistoGUI_RCS} )

# this will run uic on .ui files:
QT4_WRAP_UI( HistoGUI_UI_HDRS ${HistoGUI_UIS} )

# and finally this will run moc:
QT4_WRAP_CPP( HistoGUI_MOC_SRCS HistoGUI.h )

ADD_EXECUTABLE(HistoGUI main.cpp
             ${HistoGUI_SRCS}
             ${HistoGUI_HDRS} 
             ${HistoGUI_MOC_SRCS} 
             ${HistoGUI_UI_HDRS}
             ${HistoGUI_RC_SRCS}
)

TARGET_LINK_LIBRARIES( HistoGUI ${QT_LIBRARIES} QVTK ftkImage TinyXML libSVM vtkVolumeRendering NuclearSegmentation ftkGUI)

#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
#SET_TARGET_PROPERTIES( HistoGUI PROPERTIES WIN32_EXECUTABLE 1)
#OR MAKE A CONSOLE APPLICATION
SET_TARGET_PROPERTIES( HistoGUI PROPERTIES WIN32_EXECUTABLE 0)

INSTALL( TARGETS HistoGUI DESTINATION bin )


