# set project's name
PROJECT( Rendering )

ADD_EXECUTABLE(render main.cpp)
TARGET_LINK_LIBRARIES(render vtkRendering vtkCommon vtkIO vtksys vtkFiltering vtkHybrid)
TARGET_LINK_LIBRARIES(render      ${ITK_LIBRARIES} )
TARGET_LINK_LIBRARIES(render tiff Trace TinyXML ${QT_LIBRARIES})
INSTALL( TARGETS render DESTINATION ${INSTALL_BIN_DIR} )

