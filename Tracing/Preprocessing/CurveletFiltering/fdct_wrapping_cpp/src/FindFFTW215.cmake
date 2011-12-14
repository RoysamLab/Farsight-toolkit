
SET(FFTW_INCLUDE_SEARCH_PATH
	/usr/include
	/usr/include/fftw
	/usr/local/include
	/usr/local/include/fftw
	C:/users/Arun/Research/fftw/fftw_2.1.5_install/include
	${FFTW215_INSTALL_DIR}/include
	)

SET(FFTW_LIB_SEARCH_PATH
	/usr/lib
	/usr/lib/fftw
	/usr/local/lib
	/usr/local/lib/fftw
	C:/users/Arun/Research/fftw/fftw_2.1.5_install/lib
	${FFTW215_INSTALL_DIR}/lib
	)

FIND_PATH(FFTW_INCLUDE_PATH fftw.h ${FFTW_INCLUDE_SEARCH_PATH})

IF(FFTW_INCLUDE_PATH)
	SET(FFTW_INCLUDE_DIR ${FFTW_INCLUDE_PATH})
ENDIF(FFTW_INCLUDE_PATH)


FIND_LIBRARY(FFTW_LIB_SEARCH 
	NAMES libfftw215 libfftw fftw 
	PATHS ${FFTW_LIB_SEARCH_PATH})


FIND_LIBRARY(RFFTW_LIB_SEARCH 
	NAMES librfftw215 libfftw fftw
	PATHS ${FFTW_LIB_SEARCH_PATH})

FIND_LIBRARY(FFTW_THREADS_LIB_SEARCH 
	NAMES libfftw215_threads fftw_threads libfftw_threads
	PATHS ${FFTW_LIB_SEARCH_PATH})


IF(FFTW_LIB_SEARCH) #assuming that if you find one, you'd find the rest in the same place
	SET(FFTW_LIBRARIES ${FFTW_LIB_SEARCH} ${RFFTW_LIB_SEARCH} ${FFTW_THREADS_LIB_SEARCH})
	SET(FFTW_FOUND TRUE)
	GET_FILENAME_COMPONENT(FFTW_LINK_DIRECTORIES ${FFTW_LIB_SEARCH} PATH)
	#MESSAGE("Debug: FFTW_LIBRARIES=${FFTW_LIBRARIES}")
	#MESSAGE("Debug: FFTW_LINK_DIRECTORIES=${FFTW_LINK_DIRECTORIES}")
ELSE(FFTW_LIB_SEARCH)
	MESSAGE("Could not find FFTW anywhere.")
ENDIF(FFTW_LIB_SEARCH)


MARK_AS_ADVANCED(
	FFTW_INCLUDE_DIR
	FFTW_LIBRARIES
	FFTW_LINK_DIRECTORIES)

