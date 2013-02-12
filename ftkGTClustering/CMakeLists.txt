cmake_minimum_required(VERSION 2.8)

PROJECT (ftkGTClustering)

INCLUDE_DIRECTORIES(${VXL_MUL_INCLUDE_DIR}) 

#FIND_PACKAGE(VXL)

SET(ftkGTClustering_sources
	ftkGTClustering.cxx
	)
	
SET(ftkGTClustering_headers
	ftkGTClustering.h
	)	

#ADD_EXECUTABLE ( ftkGTClustering ${ftkGTClustering_sources} ${ftkGTClustering_headers})
#TARGET_LINK_LIBRARIES( ftkGTClustering vcl vnl vnl_algo vnl_io mbl)

ADD_LIBRARY(ftkGTClustering_lib ${ftkGTClustering_sources} ${ftkGTClustering_headers})
TARGET_LINK_LIBRARIES(ftkGTClustering_lib ${VXL_LIBRARIES} vcl vnl vnl_algo vnl_io mbl)