PROJECT(Daniela)

#Build the Daniela class files
ADD_LIBRARY(Daniela Daniela.cpp  Daniela.h itkRGBToHSVColorSpacePixelAccessor.h)
ADD_EXECUTABLE(Daniela_main  Daniela_main.cpp  )
TARGET_LINK_LIBRARIES(Daniela_main Daniela ITKIO ITKCommon ITKBasicFilters ITKAlgorithms)

#QT STUFF
QT4_WRAP_CPP( Daniela_MOC_SRCS Daniela_GUI.h)

#GUI stuff
ADD_EXECUTABLE(Daniela_GUI main.cpp Daniela_GUI.cpp Daniela_GUI.h Daniela.h ${Daniela_MOC_SRCS})
TARGET_LINK_LIBRARIES(Daniela_GUI Daniela ${QT_LIBRARIES} ITKCommon ftkGUI)

