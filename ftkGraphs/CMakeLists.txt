SET( FTKGRAPH_SRCS
  FTKgraph.cpp
  kNearestObjects.txx
)

SET( FTKGRAPH_HDRS
  FTKgraph.h
  kNearestObjects.h
)

ADD_LIBRARY(ftkGraphs ${FTKGRAPH_SRCS} ${FTKGRAPH_HDRS})

TARGET_LINK_LIBRARIES( ftkGraphs vtkCommon vtkInfovis vtkViews vtkRendering vtkFiltering vtkGraphics   ftkGUI ftkCommon ftkFeatures     ${ITK_LIBRARIES} )
