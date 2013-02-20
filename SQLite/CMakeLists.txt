SET( SQLITE3_SRCS
  sqlite3.h
  sqlite3ext.h
)

SET( SQLITE3_HDRS
  shell.c
  sqlite3.c
  sqlite3.def
)

ADD_LIBRARY( sqlite3
            ${SQLITE3_SRCS} 
            ${SQLITE3_HDRS} 
)