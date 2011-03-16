include(${QT_USE_FILE})
INCLUDE(${ITK_USE_FILE})
INCLUDE(${VTK_USE_FILE})

SET( TracingSysLibraries
    ITKCommon
    ITKNumerics
    ITKStatistics
    ITKIO
    ITKBasicFilters
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
      New_ObjectSelection.cpp
      OpenSnakeTracer.cpp

      Rendering/New_ImageActors.h
      Rendering/New_ImageActors.cxx
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
  
      TracingCore/Vesselness/itkHessianToObjectnessMeasureImageFilter.h
      TracingCore/Vesselness/itkHessianToObjectnessMeasureImageFilter.txx
      TracingCore/Vesselness/itkMultiScaleHessianBasedMeasureImageFilter.h
      TracingCore/Vesselness/itkMultiScaleHessianBasedMeasureImageFilter.txx

      TracingCore/Vesselness/itkVesselEnhancingDiffusion3DImageFilter.h
      TracingCore/Vesselness/itkVesselEnhancingDiffusion3DImageFilter.txx
      TracingCore/Vesselness/itkAnisotropicDiffusionVesselEnhancementImageFilter.h
      TracingCore/Vesselness/itkAnisotropicDiffusionVesselEnhancementImageFilter.txx
      TracingCore/Vesselness/itkAnisotropicDiffusionVesselEnhancementFunction.h
      TracingCore/Vesselness/itkAnisotropicDiffusionVesselEnhancementFunction.txx
      TracingCore/Vesselness/itkHessianSmoothed3DToVesselnessMeasureImageFilter.h
      TracingCore/Vesselness/itkHessianSmoothed3DToVesselnessMeasureImageFilter.txx
      TracingCore/Vesselness/itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.h
      TracingCore/Vesselness/itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.txx
      TracingCore/Vesselness/itkSymmetricEigenVectorAnalysisImageFilter.h

      TracingCore/GCv2p2/GCoptimization.cpp
      TracingCore/GCv2p2/GCoptimization.h
      TracingCore/GCv2p2/graph.cpp
      TracingCore/GCv2p2/graph.h
      TracingCore/GCv2p2/LinkedBlockList.cpp
      TracingCore/GCv2p2/maxflow.cpp

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
