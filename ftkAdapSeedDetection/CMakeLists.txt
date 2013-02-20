CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE( ${ITK_USE_FILE} )
ELSE(ITK_FOUND)
      MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without VTK")
ENDIF(VTK_FOUND)

LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS})
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS})


################################################################################
################################################################################


SET( FTKVOTING_SRCS

main.cpp
ftkVoting.cpp
ftkVoting_3D.cpp
ftkVotingGlobal.cpp
 
)

################################################################################
################################################################################


SET( FTKVOTING_HDRS

ftkVoting.h
ftkVoting_3D.h
ftkVotingGlobal.h
 
)

################################################################################
################################################################################


ADD_LIBRARY(LIBftkVoting ${FTKVOTING_HDRS} ${FTKVOTING_SRCS})

ADD_EXECUTABLE( ftkVoting
		main.cpp
)

#ADD_EXECUTABLE( ftkVoting
#		main.cpp
#		${FTKVOTING_HDRS}
#		${FTKVOTING_SRCS}
#)

################################################################################
################################################################################


TARGET_LINK_LIBRARIES ( ftkVoting ${VTK_LIBRARIES} ${ITK_LIBRARIES} LIBftkVoting )
#TARGET_LINK_LIBRARIES ( ftkVoting ${VTK_LIBRARIES} ${ITK_LIBRARIES} )
