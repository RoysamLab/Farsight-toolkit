PROJECT (ImageManager)

SET( Libraries
    ${ITK_LIBRARIES}
    #ftkCommon
    ${QT_LIBRARIES}
    ftkGUI )
	
	QT4_WRAP_CPP( GUI_MOC_SRCS ImageManager.h)
	
ADD_LIBRARY( ImageManagerSources
	ImageManager.h	ImageManager.cpp
  ${GUI_MOC_SRCS} 
	
)

TARGET_LINK_LIBRARIES( ImageManagerSources ${Libraries})
ADD_EXECUTABLE(${PROJECT_NAME} ${EXE_TYPE} ImageManagerDriver.cpp )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}  ImageManagerSources)

if(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION bin )
endif(APPLE)