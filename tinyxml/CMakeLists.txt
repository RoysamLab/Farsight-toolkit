SET( TINYXML_SRCS    
	tinystr.cpp
	tinyxml.cpp
	tinyxmlerror.cpp
	tinyxmlparser.cpp

)

SET( TINYXML_HDRS
	tinystr.h
	tinyxml.h
)

ADD_LIBRARY( TinyXML
			${TINYXML_SRCS}
			${TINYXML_HDRS}
)

#ADD_EXECUTABLE( XMLtoTXT XMLtoTXT.cpp)
#TARGET_LINK_LIBRARIES( XMLtoTXT TinyXML )

#ADD_EXECUTABLE( XMLtoMETA XMLtoMETA.cpp)
#TARGET_LINK_LIBRARIES( XMLtoMETA TinyXML )


#INSTALL( TARGETS TinyXML DESTINATION . )

