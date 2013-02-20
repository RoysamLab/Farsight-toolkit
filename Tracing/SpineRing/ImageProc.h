#ifndef IMAGEPROC_H
#define IMAGEPROC_H
#include "SpineRing.h"


// templating these functions was a pain bec. of compiler 
// warnings on explicit instantiations. Decided to write
// them explicitly overloaded for now. Maybe a better compiler 
// later on will not be as ...
SpineImage2DType::Pointer MIPGenerator(SpineImageType::Pointer im);
SpeedImage2DType::Pointer MaxIPGenerator(SpeedImageType::Pointer im);
SpineImage2DType::Pointer MaxIPGenerator(SpineImageType::Pointer im);
int writeImageToFile (SpineImage2DType::Pointer,  std::string);
int writeImageToFile (SpeedImage2DType::Pointer,  std::string);
int writeImageToFile (SpeedImageType::Pointer,  std::string);
int writeImageToFile (RGBImageType::Pointer,      std::string);
int writeImageToFile (BinaryImage2DType::Pointer, std::string);
int writeImageToFile (SpineImageType::Pointer, std::string);
void ImStatsSpeedImageType(SpeedImageType::Pointer im, std::string s);
//template <typename T> void ImStats(typename T::Pointer im);

//template <typename T> int writeImageToFile(typename T::Pointer im, std::string fn);

//template int writeImageToFile <SpineImage2DType> (SpineImage2DType::Pointer, std::string);

//template <typename T, typename U>
//typename U::Pointer MIPGenerator(typename T::Pointer im);

//template SpineImage2DType::Pointer MIPGenerator<SpineImageType, SpineImage2DType>(SpineImageType::Pointer im);
//template SpeedImage2DType::Pointer MIPGenerator<SpeedImageType, SpeedImage2DType>(SpeedImageType::Pointer im);


class ImageDebugger {
public:
	ImageDebugger(SpineImageType::Pointer im, std::string fn): im3D(im), MIP(0), basefilename(fn)
	{
		MIP     = MIPGenerator(im);
		rgbim2D = CreateRGBImFrom2D(MIP);
	};
	~ImageDebugger(){};
	SpineImageType::Pointer im3D;
	SpineImage2DType::Pointer MIP;
	//BinaryImageType::Pointer bin3D;
	BinaryImage2DType::Pointer bin2D;
	RGBImageType::Pointer rgbim2D;

    void ImageStatistics(SpineImageType::Pointer & im3D);
	//void MIPGenerator();
	void MaxIPGenerator(BinaryImageType::Pointer bin3D);
	void WriteTempObj2D(int trnum, int ringnum, IndexVecType	*inring_f_idx);
	void WriteSpCandidates2RGBImage(GlobalDetector *gd);
	//void Write2DRGBImageToFile(RGBImageType::Pointer rgbim2D, std::string fn);
	
	RGBImageType::Pointer CreateRGBImFrom2D(SpineImage2DType::Pointer im2d);


	int  PlotContOnRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color, PointSetContainerType::Pointer pcont);
	//PointSetContainerType::Pointer MakeBox(PointType *pt);
	void PlotPtAsBoxOnRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color, RGBPointType pt, int plotboxx, int plotboxy);
	void PlotPtAsBoxOnRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color, PointType    pt, int plotboxx, int plotboxy);
	void PlotContAsBoxOnRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color, PointSetContainerType::Pointer pcont, int plotboxx, int plotboxy);

	std::string basefilename;
	bool NeedMIP;
};

#endif
