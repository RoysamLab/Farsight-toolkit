ADD_EXECUTABLE ( segren segren.cpp )
TARGET_LINK_LIBRARIES ( segren ftkCommon NuclearSegmentation    vtkRendering vtkIO vtkHybrid vtkCommon vtkFiltering vtkVolumeRendering vtksys      ${ITK_LIBRARIES} )
INSTALL ( TARGETS segren DESTINATION bin )
