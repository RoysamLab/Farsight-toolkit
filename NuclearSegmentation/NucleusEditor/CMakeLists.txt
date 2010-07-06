SET( NucleusEditor_SRCS
  NucleusEditor.cpp
  ProjectFilenamesDialog.cpp
  ExclusionDialog.cpp
  PreferencesDialog.cpp
  ftkProjectDefinition.cpp
  ftkProjectProcessor.cpp
  ftkProjectFiles.cpp
)

SET( NucleusEditor_HDRS
  NucleusEditor.h
  ProjectFilenamesDialog.h
  ExclusionDialog.h
  PreferencesDialog.h
  ftkProjectDefinition.h
  ftkProjectProcessor.h
  ftkProjectFiles.h
)

SET( NucleusEditor_MOC_REQD_HDRS
  NucleusEditor.h 
  ProjectFilenamesDialog.h 
  ExclusionDialog.h
  PreferencesDialog.h
)

SET( NucleusEditor_UIS
)

SET( NucleusEditor_RCS
  icons.qrc
)

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

TARGET_LINK_LIBRARIES( NucleusEditor ${QT_LIBRARIES} ftkGUI ftkImage ftkCommon ftkFeatures libSVM vtkVolumeRendering QVTK NuclearSegmentation CytoplasmSegmentation PixelLevelAnalysis )

ADD_EXECUTABLE( projproc processor_main.cpp ftkProjectProcessor.cpp ftkProjectDefinition.cpp)
TARGET_LINK_LIBRARIES( projproc ftkImage ftkCommon ftkFeatures NuclearSegmentation CytoplasmSegmentation PixelLevelAnalysis )
install(TARGETS projproc DESTINATION bin )

#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
#SET_TARGET_PROPERTIES( NucleusEditor PROPERTIES WIN32_EXECUTABLE 1)
#OR MAKE A CONSOLE APPLICATION
SET_TARGET_PROPERTIES( NucleusEditor PROPERTIES WIN32_EXECUTABLE 0)

if(APPLE)
   install(TARGETS NucleusEditor DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS NucleusEditor DESTINATION bin )
endif(APPLE)

