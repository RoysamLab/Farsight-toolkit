CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


PROJECT(ExtractSomas)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

INCLUDE(CTest)

ADD_EXECUTABLE(ImageHistogram2  ImageHistogram2.cxx)
ADD_EXECUTABLE(KappaSigmaThresholdImageFilterTest KappaSigmaThresholdImageFilterTest.cxx)
ADD_EXECUTABLE(RobustAutomaticThresholdImageFilterTest RobustAutomaticThresholdImageFilterTest.cxx)
ADD_EXECUTABLE(RobustAutomaticThresholdImageFilterTest2D RobustAutomaticThresholdImageFilterTest2D.cxx)
ADD_EXECUTABLE(MorphologyGrayscaleFilters MorphologyGrayscaleFilters.cxx)
ADD_EXECUTABLE(ImageReadWrite ImageReadWrite.cxx)
ADD_EXECUTABLE(MedianImageFilter MedianImageFilter.cxx)
ADD_EXECUTABLE(PadImage PadImage.cxx)
ADD_EXECUTABLE(PadDistanceMap PadDistanceMap.cxx)
ADD_EXECUTABLE(HoleFillingFilterTest HoleFillingFilterTest.cxx)
ADD_EXECUTABLE(ExtractImageFromPadding ExtractImageFromPadding.cxx)
ADD_EXECUTABLE(DistanceMapFilter DistanceMapFilter.cxx)
ADD_EXECUTABLE(BinaryThresholdImageFilter BinaryThresholdImageFilter.cxx)
ADD_EXECUTABLE(DownwardFrontPropagation DownwardFrontPropagation.cxx)
ADD_EXECUTABLE(VoronoiFilter VoronoiFilter.cxx)

TARGET_LINK_LIBRARIES(ImageHistogram2 ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(KappaSigmaThresholdImageFilterTest ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(RobustAutomaticThresholdImageFilterTest ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(RobustAutomaticThresholdImageFilterTest2D ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(MorphologyGrayscaleFilters ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(ImageReadWrite ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(MedianImageFilter ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(PadImage ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(PadDistanceMap ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(HoleFillingFilterTest ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(ExtractImageFromPadding ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(DistanceMapFilter ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(BinaryThresholdImageFilter ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(DownwardFrontPropagation ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(VoronoiFilter ITKIO ITKStatistics)

ADD_TEST(Read ImageReadWrite /path/to/image.tif input.tif )
ADD_TEST(Median MedianImageFilter input.tif median.tif)
    SET_TESTS_PROPERTIES ( Median PROPERTIES DEPENDS Read )
ADD_TEST(Robust RobustAutomaticThresholdImageFilterTest median.tif robust.tif 2 10.0)
    SET_TESTS_PROPERTIES ( Robust PROPERTIES DEPENDS Median )
ADD_TEST(Pad1 PadImage robust.tif padded.tif)
    SET_TESTS_PROPERTIES ( Pad1 PROPERTIES DEPENDS Robust )
ADD_TEST(HoleFilling HoleFillingFilterTest padded.tif filled_padded.tif 1000 1 )
    SET_TESTS_PROPERTIES ( HoleFilling PROPERTIES DEPENDS Pad1 )
ADD_TEST(Depad1 ExtractImageFromPadding filled_padded.tif filled.tif)
    SET_TESTS_PROPERTIES ( Depad1 PROPERTIES DEPENDS HoleFilling )
ADD_TEST(Distance DistanceMapFilter filled.tif distance.mhd)
    SET_TESTS_PROPERTIES ( Distance PROPERTIES DEPENDS Depad1 )

ADD_TEST(Threshold BinaryThresholdImageFilter distance.mhd distanceThreshold.tif 
  5.0  # 80% of expected blob size
  )
    SET_TESTS_PROPERTIES ( Threshold PROPERTIES DEPENDS Distance )

ADD_TEST(Pad2 PadImage distanceThreshold.tif distanceThreshold_padded.tif)
    SET_TESTS_PROPERTIES ( Pad2 PROPERTIES DEPENDS Threshold )

ADD_TEST(PadDistanceMap PadDistanceMap distance.mhd distance_padded.mhd)
    SET_TESTS_PROPERTIES ( PadDistanceMap PROPERTIES DEPENDS Pad2 )

ADD_TEST(Propagation DownwardFrontPropagation 
  distanceThreshold_padded.tif 
  distance_padded.mhd
  output_somas_padded.tif
  1000
  4.0  # expected radius of tubes
  )
    SET_TESTS_PROPERTIES ( Propagation PROPERTIES DEPENDS PadDistanceMap )

ADD_TEST(Depad2 ExtractImageFromPadding output_somas_padded.tif output_somas.tif)
    SET_TESTS_PROPERTIES ( Depad2 PROPERTIES DEPENDS Propagation )

#ADD_TEST(Voronoi1 VoronoiFilter somas.tif voronoi.mhd)
#    SET_TESTS_PROPERTIES ( Voronoi1 PROPERTIES DEPENDS Propagation )
