SET( NESQLITE_SRCS    
	NESqliteFactory.cpp
)

SET( NESQLITE_HDRS
	NESqliteFactory.h
)

ADD_LIBRARY( NESQLITE
             ${NESQLITE_SRCS} 
             ${NESQLITE_HDRS} 
)

TARGET_LINK_LIBRARIES( sqlite3 )

INSTALL(FILES NE.s3db DESTINATION bin )