# set project's name
PROJECT( Rendering )

ADD_EXECUTABLE(render main.cpp)
TARGET_LINK_LIBRARIES(render vtkRendering vtkCommon vtkIO vtksys vtkFiltering vtkHybrid)
TARGET_LINK_LIBRARIES(render ITKBasicFilters ITKCommon ITKIO)
TARGET_LINK_LIBRARIES(render Trace TinyXML)
INSTALL(TARGETS render DESTINATION bin)
