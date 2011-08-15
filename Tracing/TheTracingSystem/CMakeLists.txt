SET( TracingSysLibraries
    ${ITK_LIBRARIES}
    rtvl
    rgtl
    vtkRendering
    vtkIO
    vtkHybrid
    vtkWidgets
    vtkCommon
    vtksys
    vtkFiltering
    vtkVolumeRendering
    QVTK
)

SET( GUI_RCS
  icons.qrc
)

QT4_ADD_RESOURCES( GUI_RC_SRCS ${GUI_RCS} )
QT4_WRAP_CPP(GUI_MOC_SRCS TracingGUI.h TracingView.h ObjectSelection.h ParametersGroup.h SlidersGroup.h OpenSnakeTracer.h dialogs.h dialogs_montage.h)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )


ADD_EXECUTABLE(TracingSystem 
      main.cpp 
      SlidersGroup.cpp
      ParametersGroup.cpp
      dialogs.cpp 
      dialogs_montage.cpp
      TracingGUI.cpp
      TracingView.cpp  
      ObjectSelection.cpp
      OpenSnakeTracer.cpp

      Rendering/ImageActors.h
      Rendering/ImageActors.cxx
      Rendering/itkImageToVTKImageFilter.h
      Rendering/itkImageToVTKImageFilter.txx

      TracingCore/PointOperation.h
      TracingCore/PointOperation.cpp
      TracingCore/ImageOperation.h
      TracingCore/ImageOperation.cpp
      TracingCore/SnakeOperation.h
      TracingCore/SnakeOperation.cpp

      TracingCore/Filters/itkLabelGeometryImageFilter.h
      TracingCore/Filters/itkLabelGeometryImageFilter.txx
      TracingCore/Filters/itkBinaryThinningImageFilter3D.h
      TracingCore/Filters/itkBinaryThinningImageFilter3D.txx
      TracingCore/Filters/itkMaximumProjectionImageFilter.h
      TracingCore/Filters/itkProjectionImageFilter.h
      TracingCore/Filters/itkProjectionImageFilter.txx
  
      TracingCore/Montage/Register.h
      TracingCore/Montage/Register.cpp

     ${GUI_MOC_SRCS}
     ${GUI_RC_SRCS})


TARGET_LINK_LIBRARIES(TracingSystem ${QT_LIBRARIES} ${TracingSysLibraries} TinyXML)
 
#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
INSTALL( TARGETS TracingSystem DESTINATION bin )

#INCLUDE(InstallRequiredSystemLibraries)
if(APPLE)
   install(TARGETS TracingSystem DESTINATION ${FARSIGHT_BUNDLE_LOCATION} )
else(APPLE)
   install(TARGETS TracingSystem DESTINATION bin )
endif(APPLE)

