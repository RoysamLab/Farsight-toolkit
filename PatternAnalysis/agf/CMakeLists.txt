SET( AGF_SRCS
	#badri.cpp
	agf.cpp
)

SET( AGF_HDRS
	agf.h
)


ADD_LIBRARY( libagf
            ${AGF_SRCS}
            ${AGF_HDRS} 
)

