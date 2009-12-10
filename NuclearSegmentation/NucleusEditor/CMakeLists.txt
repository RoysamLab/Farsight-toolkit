SET( NucluesEditor_SRCS
  NucleusEditor.cpp
  ProjectFilenamesDialog.cpp
  ftkProjectDefinition.cpp
  ftkProjectProcessor.cpp
  ftkProjectFiles.cpp
)

SET( NucluesEditor_HDRS
  NucleusEditor.h
  ProjectFilenamesDialog.h
  ftkProjectDefinition.h
  ftkProjectProcessor.h
  ftkProjectFiles.h
)

SET( NucluesEditor_UIS
)

SET( NucluesEditor_RCS
  icons.qrc
)

# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
QT4_ADD_RESOURCES( NucluesEditor_RC_SRCS ${NucluesEditor_RCS} )

# this will run uic on .ui files:
QT4_WRAP_UI( NucluesEditor_UI_HDRS ${NucluesEditor_UIS} )

# and finally this will run moc:
QT4_WRAP_CPP( NucluesEditor_MOC_SRCS NucleusEditor.h ProjectFilenamesDialog.h )

ADD_EXECUTABLE(NucleusEditor main.cpp
             ${NucluesEditor_SRCS}
             ${NucluesEditor_HDRS} 
             ${NucluesEditor_MOC_SRCS} 
             ${NucluesEditor_UI_HDRS}
             ${NucluesEditor_RC_SRCS}
)

TARGET_LINK_LIBRARIES( NucleusEditor ${QT_LIBRARIES} ftkGUI ftkImage ftkCommon libSVM vtkVolumeRendering QVTK NuclearSegmentation CytoplasmSegmentation )

#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
#SET_TARGET_PROPERTIES( NucleusEditor PROPERTIES WIN32_EXECUTABLE 1)
#OR MAKE A CONSOLE APPLICATION
SET_TARGET_PROPERTIES( NucleusEditor PROPERTIES WIN32_EXECUTABLE 0)

INSTALL( TARGETS NucleusEditor DESTINATION bin )


