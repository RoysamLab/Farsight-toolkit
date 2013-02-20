PROJECT(Dendrogram)
SET( Dendo_SRCS
  
  Dendrogram.cpp
  ClusClus.cpp
)

SET( Dendo_HDRS
  
  Dendrogram.h
  ClusClus.h
)

SET( Dendo_MOC_REQD_HDRS
	Dendrogram.h

	)

QT4_WRAP_CPP( Dendo_MOC_REQD_HDRS Dendro_main.cpp Dendrogram.h)


ADD_EXECUTABLE( Dendrogram Dendro_main.cpp
               
                ${Dendo_SRCS}
				${Dendo_HDRS}
				${Dendo_MOC_REQD_HDRS}
                
                )
                
ADD_LIBRARY( Dendrogram_lib ${Dendo_SRCS} ${Dendo_HDRS} ${Dendo_MOC_REQD_HDRS})

TARGET_LINK_LIBRARIES( Dendrogram
  ${VTK_LIBRARIES}
  ${QT_LIBRARIES}
  ftkGUI
  )
  
TARGET_LINK_LIBRARIES( Dendrogram_lib ${QT_LIBRARIES} )

# and finally this will run moc:
