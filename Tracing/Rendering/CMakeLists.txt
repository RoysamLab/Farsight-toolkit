# set project's name
PROJECT( Rendering )

ADD_EXECUTABLE(render main.cpp)
TARGET_LINK_LIBRARIES(render
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  tiff
  Trace
  TinyXML
  ${QT_LIBRARIES}
  )
INSTALL( TARGETS render DESTINATION ${INSTALL_BIN_DIR} )

