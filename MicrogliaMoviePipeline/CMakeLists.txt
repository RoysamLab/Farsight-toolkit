cmake_minimum_required(VERSION 2.8)
PROJECT(MicrogliaMoviePipeline)
ADD_EXECUTABLE(MicrogliaMoviePipeline movie_main.cpp SegAndTrace.h SegAndTrace.cpp)
TARGET_LINK_LIBRARIES(MicrogliaMoviePipeline ${ITK_LIBRARIES} Yousef_Nucleus_Seg segment_soma MultipleNeuronTracerLib)
