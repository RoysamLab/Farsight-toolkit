INCLUDE_DIRECTORIES("${Farsight_SOURCE_DIR}/ftkCommon")
INCLUDE_DIRECTORIES( ${Farsight_BINARY_DIR}/Tracing/MDL2/wizard ${Farsight_SOURCE_DIR}/Tracing/MDL2/wizard )

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
                      mdlIntegratedSkeleton mdlMST mdlBSplineFitting TinyXML)

