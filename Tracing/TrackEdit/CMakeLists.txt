CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

# set project's name
PROJECT( TrackingEditor )

 

SET(Libraries
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
)

QT4_WRAP_CPP( GUI_MOC_SRCS TrackingImageView.h TrackingEventHandler.h TrackSelection.h TrackingModel.h TrackingVolumeView.h TrackingKymoView.h)

  ADD_EXECUTABLE(TrackingEditor main.cpp TrackSelection.cpp TrackSelection.h 
	  TrackingImageView.cpp TrackingImageView.h TrackingEventHandler.h TrackingEventHandler.cpp TrackingModel.h 
	  TrackingModel.cpp TrackingVolumeView.cpp TrackingVolumeView.h TrackingKymoView.h TrackingKymoView.cpp Trace.h Trace.cpp ${GUI_MOC_SRCS})
  ADD_LIBRARY(helpers helpers.cpp helpers.h)
  TARGET_LINK_LIBRARIES(TrackingEditor ${Libraries} ${QT_LIBRARIES} helpers)
