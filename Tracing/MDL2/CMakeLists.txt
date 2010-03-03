PROJECT(MDL2)
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} ${Farsight_BINARY_DIR}/Tracing/MDL2 ${Farsight_SOURCE_DIR}/Tracing/MDL2 )

SET( cbin_SRCS
  GCBinarization/cell_binarization.h
  GCBinarization/cell_binarization.cpp    
  GCBinarization/new_graph.h
  GCBinarization/block_b.h   
)

#CHECK TO MAKE SURE ITK LIBRARIES ARE NEEDED FOR THESE:

ADD_LIBRARY(mdlVolumeProcess mdlVolumeProcess.h mdlVolumeProcess.cpp ${cbin_SRCS})
TARGET_LINK_LIBRARIES(mdlVolumeProcess ITKIO ITKStatistics ITKBasicFilters ITKCommon)

ADD_LIBRARY(mdlIntegratedSkeleton mdlIntegratedSkeleton.h mdlIntegratedSkeleton.cpp)
TARGET_LINK_LIBRARIES(mdlVolumeProcess ITKBasicFilters ITKCommon)

ADD_LIBRARY(mdlUtils mdlUtils.h mdlUtils.cpp)

ADD_LIBRARY(mdlMST mdlMST.h mdlMST.cpp WeightedMahalsnobisDistance.cxx)
TARGET_LINK_LIBRARIES(mdlMST mdlUtils)

ADD_LIBRARY(mdlBSplineFitting mdlBSplineFitting.cpp mdlBSplineFitting.h)

ADD_EXECUTABLE(volumeProcess2 main.cpp)
TARGET_LINK_LIBRARIES(volumeProcess2 mdlVolumeProcess)

ADD_EXECUTABLE(mdlpipe pipeline_main.cpp)
TARGET_LINK_LIBRARIES(mdlpipe mdlVolumeProcess mdlIntegratedSkeleton mdlMST mdlBSplineFitting  )

QT4_WRAP_UI( Wizard_HDRS MDL2Wizard.ui )
QT4_WRAP_CPP( Wizard_MOC_SRCS
  MDL2Wizard.h MDL2WizardHelper.h
  MDL2IntroPage.h MDL2PreprocessingPage.h
  MDL2PhaseOnePage.h MDL2PhaseTwoPage.h
)
QT4_ADD_RESOURCES(
  Wizard_RCC_SRCS
  MDL2Wizard.qrc
)
SET(MDL2Wizard_SRCS
  MDL2Wizard.h             MDL2Wizard.cxx 
  MDL2WizardHelper.h       MDL2WizardHelper.cxx 
  MDL2IntroPage.h          MDL2IntroPage.cxx
  MDL2PreprocessingPage.h  MDL2PreprocessingPage.cxx
  MDL2PhaseOnePage.h       MDL2PhaseOnePage.cxx
  MDL2PhaseTwoPage.h       MDL2PhaseTwoPage.cxx
  MDL2Wizard-main.cxx
  ${Wizard_MOC_SRCS} ${Wizard_HDRS} ${Wizard_RCC_SRCS}
)
ADD_EXECUTABLE(MDL2Wizard ${EXE_TYPE} ${MDL2Wizard_SRCS})
TARGET_LINK_LIBRARIES(MDL2Wizard ${QT_LIBRARIES} vtkCommon vtkIO vtkRendering
                      vtkVolumeRendering QVTK mdlVolumeProcess
                      mdlIntegratedSkeleton mdlMST mdlBSplineFitting  )

