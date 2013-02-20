SET( MontagingToTraceEditorXMLConverter_SRCS
	MontagingToTraceEditorXMLConverter.cpp
	main.cpp
)

SET( MontagingToTraceEditorXMLConverter_HDRS
	MontagingToTraceEditorXMLConverter.h
	ImageTransformData.h
)

#MESSAGE(STATUS "Creating MontagingToTraceEditorXMLConverter Project")

#LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS})
#INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS})
#TARGET_LINK_LIBRARIES( )

ADD_EXECUTABLE( MontagingToTraceEditorXMLConverter ${MontagingToTraceEditorXMLConverter_HDRS} ${MontagingToTraceEditorXMLConverter_SRCS} )
TARGET_LINK_LIBRARIES( MontagingToTraceEditorXMLConverter TinyXML )

INSTALL( TARGETS MontagingToTraceEditorXMLConverter RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
