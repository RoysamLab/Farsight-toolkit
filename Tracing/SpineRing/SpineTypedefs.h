#ifndef SPINETYPEDEFS_H
#define SPINETYPEDEFS_H


// No need to change these for FTK
typedef itk::Image< unsigned short, spr_SPIMGDIMENSION>   CCOutputImageType;

//FTK FIXME: this needs to be modified when part of FTK
typedef itk::ImageFileReader< SpineImageType > ReaderType;

typedef itk::ConnectedComponentImageFilter< BinaryImageType, CCOutputImageType > CCFilterType;
typedef itk::RelabelComponentImageFilter< CCOutputImageType, CCOutputImageType > RelabelType;

//typedef itk::LinearInterpolateImageFunction<SpineImageType, double>  InterpolatorType;
//typedef itk::NearestNeighborInterpolateImageFunction<ImageType3D,float>  InterpolatorType;
//typedef InterpolatorType::PointType InterpPointType;


typedef itk::RGBPixel<unsigned char>   RGBPixelType;
typedef itk::Image< RGBPixelType,  2 >    RGBImageType;
typedef RGBImageType::IndexType RGBIndexType;
typedef	RGBImageType::PointType RGBPointType;

typedef float											SpeedPixelType; 
typedef itk::Image< SpeedPixelType, spr_SPIMGDIMENSION> SpeedImageType;
typedef itk::Image< SpeedPixelType, 2> SpeedImage2DType;
typedef itk::PolyLineParametricPath< spr_SPIMGDIMENSION > PathType;

class SpCandidate;
typedef std::map<unsigned short, SpCandidate*> SpCandMapType;

#endif
