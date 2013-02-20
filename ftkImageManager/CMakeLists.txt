PROJECT (ImageManager)

SET( Libraries
    ${ITK_LIBRARIES}
    ${QT_LIBRARIES}
    ftkGUI )
	
QT4_WRAP_CPP( GUI_MOC_SRCS ImageManager.h)
	
ADD_LIBRARY( ImageManagerSources
	ImageManager.h	ImageManager.cpp
  ${GUI_MOC_SRCS} 
)
TARGET_LINK_LIBRARIES( ImageManagerSources ${Libraries})

add_farsight_executable(ImageManager ImageManagerDriver.cpp ImageManagerSources)

