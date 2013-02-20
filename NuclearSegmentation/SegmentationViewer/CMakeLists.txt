ADD_EXECUTABLE ( segren segren.cpp )
TARGET_LINK_LIBRARIES ( segren
  ftkCommon
  NuclearSegmentation
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  )
INSTALL ( TARGETS segren DESTINATION ${INSTALL_BIN_DIR} )
