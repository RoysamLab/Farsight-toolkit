PROJECT(MDL)
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} ${Farsight_BINARY_DIR}/Tracing/MDL 
                     ${Farsight_SOURCE_DIR}/Tracing/MDL )
SET( cbin_SRCS
  GCBinarization/cell_binarization.h
  GCBinarization/cell_binarization.cpp    
  GCBinarization/new_graph.h
  GCBinarization/block_b.h   
)
ADD_EXECUTABLE(AnisoDiffuse AnisoDiffuse.cpp)
ADD_EXECUTABLE(ConnCompntwFldfill ConnCompntwFldfill.cpp)
ADD_EXECUTABLE(GradientVecField GradientVecField.cpp)
ADD_EXECUTABLE(MinSpanTree distTransform.cpp MinSpanTree.cpp morphGraphPrune.cpp) 
ADD_EXECUTABLE(Integratedskel Integratedskel.cpp)
ADD_EXECUTABLE(ObjectVolumeProcess ObjectVolumeProcess.cpp)
TARGET_LINK_LIBRARIES(ObjectVolumeProcess       ${ITK_LIBRARIES} )
ADD_EXECUTABLE(volumeProcess distTransform.cpp volumeProcess.cpp ${cbin_SRCS})
TARGET_LINK_LIBRARIES(volumeProcess         ${ITK_LIBRARIES} )
ADD_EXECUTABLE(RefiningSkeleton RefiningSkeleton.cpp)
ADD_EXECUTABLE(MDABasedSpineExtraction distTransform.cpp MDABasedSpineExtraction.cpp morphGraphPrune.cpp WeightedMahalsnobisDistance.cxx)
TARGET_LINK_LIBRARIES(MDABasedSpineExtraction          ${ITK_LIBRARIES} )
ADD_EXECUTABLE(BackboneExtract BackboneExtract.cpp)
ADD_EXECUTABLE(ShapeDetectionLevelSegmentation ShapeDetectionLevelSetFilter.cxx)
TARGET_LINK_LIBRARIES(ShapeDetectionLevelSegmentation          ${ITK_LIBRARIES} )
ADD_EXECUTABLE(DistanceFieldComputing distTransform.cpp DistanceFieldComputing.cpp) 
TARGET_LINK_LIBRARIES(DistanceFieldComputing        ${ITK_LIBRARIES} )
ADD_EXECUTABLE(MSHS3DToVesseFilter itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilterTest.cxx)
TARGET_LINK_LIBRARIES(MSHS3DToVesseFilter          ${ITK_LIBRARIES} )

ADD_EXECUTABLE(BSplineFitting BSplineFitting.cpp)
TARGET_LINK_LIBRARIES(BSplineFitting
  ${ITK_LIBRARIES})

QT4_WRAP_UI( Wizard_HDRS MDLWizard.ui )
QT4_WRAP_CPP( Wizard_MOC_SRCS MDLWizard.h
  MDLIntroPage.h MDLPreprocessingPage.h
  MDLPhaseOnePage.h MDLPhaseTwoPage.h
)
QT4_ADD_RESOURCES(
  Wizard_RCC_SRCS
  MDLWizard.qrc
)
SET(MDLWizard_SRCS
  MDLWizard.h MDLWizard.cxx 
  MDLIntroPage.h MDLIntroPage.cxx
  MDLPreprocessingPage.h MDLPreprocessingPage.cxx
  MDLPhaseOnePage.h MDLPhaseOnePage.cxx
  MDLPhaseTwoPage.h MDLPhaseTwoPage.cxx
  MDLWizard-main.cxx
  ${Wizard_MOC_SRCS} ${Wizard_HDRS} ${Wizard_RCC_SRCS}
)
ADD_EXECUTABLE(MDLWizard ${EXE_TYPE} ${MDLWizard_SRCS})
TARGET_LINK_LIBRARIES(MDLWizard
  ${QT_LIBRARIES}
  ${VTK_LIBRARIES}
)

SET(CPACK_PACKAGE_EXECUTABLES "MDLWizard;MDL Wizard")

if(APPLE)
  SET(EXECUTABLE_LOCATION ${FARSIGHT_BUNDLE_LOCATION}/MDLWizard.app/Contents/bin)
  INSTALL( TARGETS MDLWizard
           DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
  INSTALL( TARGETS AnisoDiffuse
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS ConnCompntwFldfill
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS GradientVecField
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS MinSpanTree
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS Integratedskel
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS volumeProcess
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS ObjectVolumeProcess
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS BSplineFitting
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS RefiningSkeleton
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS MDABasedSpineExtraction
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS BackboneExtract
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS ShapeDetectionLevelSegmentation
           DESTINATION ${EXECUTABLE_LOCATION} )
  INSTALL( TARGETS DistanceFieldComputing
           DESTINATION ${EXECUTABLE_LOCATION} )
else(APPLE)
  INSTALL( TARGETS AnisoDiffuse DESTINATION bin )
  INSTALL( TARGETS ConnCompntwFldfill DESTINATION bin )
  INSTALL( TARGETS GradientVecField DESTINATION bin )
  INSTALL( TARGETS MinSpanTree DESTINATION bin )
  INSTALL( TARGETS Integratedskel DESTINATION bin )
  INSTALL( TARGETS volumeProcess DESTINATION bin )
  INSTALL( TARGETS ObjectVolumeProcess DESTINATION bin )
  INSTALL( TARGETS BSplineFitting DESTINATION bin )
  INSTALL( TARGETS RefiningSkeleton DESTINATION bin )
  INSTALL( TARGETS MDLWizard DESTINATION bin )
  INSTALL( TARGETS MDABasedSpineExtraction DESTINATION bin)
  INSTALL( TARGETS ShapeDetectionLevelSegmentation DESTINATION bin)
  INSTALL( TARGETS DistanceFieldComputing DESTINATION bin)
endif(APPLE)

