SET( SPH_SRCS 
	Spherical_Harmonic_transform.cpp
	)

SET( SPH_HDRS 
	Spherical_Harmonic_transform.h
	)

ADD_LIBRARY(libSPH ${SPH_SRCS} ${SPH_HDRS})
ADD_EXECUTABLE( SPH SPH_test.cpp )
TARGET_LINK_LIBRARIES(libSPH vnl vnl_algo vnl_io mbl libSPH)
  
TARGET_LINK_LIBRARIES( SPH libSPH)
