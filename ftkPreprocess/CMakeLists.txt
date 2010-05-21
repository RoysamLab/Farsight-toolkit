ADD_LIBRARY( ftkPreprocess
	ftkPreprocess2.h
	ftkPreprocess2.cpp
	itkBinaryThinningImageFilter3D.h
)
TARGET_LINK_LIBRARIES( ftkPreprocess ITKCommon ITKBasicFilters ITKIO ITKStatistics TinyXML)

ADD_EXECUTABLE( prep preprocess_main.cpp )
TARGET_LINK_LIBRARIES( prep ftkPreprocess )

