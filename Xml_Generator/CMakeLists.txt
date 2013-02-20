SET( XmlGen_SRCS
  Xml_Generator.cpp
)

SET( XmlGen_HDRS
  Xml_Generator.h
)

SET( XmlGen_MOC_REQD_HDRS
  Xml_Generator.h
)

# and finally this will run moc:
QT4_WRAP_CPP( XmlGen_MOC_SRCS ${XmlGen_MOC_REQD_HDRS} )

ADD_EXECUTABLE( XmlGenerator ${EXE_TYPE}
                xml_main.cpp
                ${XmlGen_SRCS}
                ${XmlGen_HDRS} 
                ${XmlGen_MOC_SRCS} 
                )

TARGET_LINK_LIBRARIES( XmlGenerator ${QT_LIBRARIES} )


