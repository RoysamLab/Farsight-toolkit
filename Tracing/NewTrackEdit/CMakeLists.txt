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

QT4_WRAP_CPP( GUI_MOC_SRCS TrackingDataModel.h)#TrackingEventHandler.h TrackSelection.h TrackingModel.h TrackingVolumeView.h TrackingKymoView.h)
#QT4_WRAP_CPP( GUI_MOC_SRCS TrackingImageView.h TrackingDataModel.h)#TrackingEventHandler.h TrackSelection.h TrackingModel.h TrackingVolumeView.h TrackingKymoView.h)

  ADD_EXECUTABLE(TrackingEditor main.cpp TrackingImageView.h TrackingImageView.cpp TrackingDataModel.h TrackingDataModel.cpp ${GUI_MOC_SRCS})
  ADD_LIBRARY(TrackEditHelper helpers.cpp)
  TARGET_LINK_LIBRARIES(TrackingEditor ${Libraries} ${QT_LIBRARIES} TrackEditHelper ftkCommon)
