SET(ftkGTClustering_sources
	ftkGTClustering.cxx
	)
	
SET(ftkGTClustering_headers
	ftkGTClustering.h
	)	

ADD_LIBRARY( ftkGTClustering_lib ${ftkGTClustering_sources} ${ftkGTClustering_headers} )
TARGET_LINK_LIBRARIES( ftkGTClustering_lib ${VXL_LIBRARIES} )
