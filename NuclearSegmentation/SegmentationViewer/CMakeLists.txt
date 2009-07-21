ADD_EXECUTABLE ( segren segren.cpp )
TARGET_LINK_LIBRARIES ( segren SegmentationCommon NuclearSegmentation ITKBasicFilters ITKCommon ITKIO vtkRendering vtkIO vtkHybrid vtkCommon vtkFiltering vtkVolumeRendering vtksys )
INSTALL ( TARGETS segren DESTINATION bin )