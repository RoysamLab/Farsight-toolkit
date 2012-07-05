
// #include "ftkMainDarpaGlobalInclude.h"

// TO ORDER

// typedef float rawImageType_flo;
typedef unsigned short Data_16_bit;
typedef unsigned char Data_8_bit;
typedef itk::Image<unsigned char,  3> rawImageType_8bit;
typedef itk::Image<Data_16_bit,  3> rawImageType_16bit;
// typedef itk::Image<Data_16_bit, 3> LabelType;
typedef itk::Image<unsigned int, 3> rawImageType_uint;
typedef itk::Image<float, 3> rawImageType_flo;
// typedef MultipleNeuronTracer::ImageType3D rawImageType_flo; // Float Image

typedef itk::ImageDuplicator< rawImageType_8bit >  rawDuplicatorType;
typedef itk::ImageDuplicator< rawImageType_16bit >  rawDuplicatorType_16;
typedef itk::ImageDuplicator< rawImageType_uint >  DuplicatorType_uint;
typedef itk::ImageDuplicator< rawImageType_flo >  DuplicatorType_flo;

typedef itk::RegionOfInterestImageFilter< rawImageType_8bit, rawImageType_8bit > rawROIFilterType_8bit;
typedef itk::RegionOfInterestImageFilter< rawImageType_16bit, rawImageType_16bit > rawROIFilterType_16bit;

typedef itk::RegionOfInterestImageFilter< rawImageType_flo, rawImageType_flo > ROIFilterType_flo;
typedef itk::RegionOfInterestImageFilter< rawImageType_uint, rawImageType_uint > ROIFilterType_uint;

typedef std::pair< rawImageType_16bit::Pointer, vtkSmartPointer< vtkTable > > segResultsType;

typedef itk::ImageFileReader<rawImageType_16bit> rawReaderType;
typedef itk::ImageFileReader<rawImageType_uint> labelReaderType;

typedef itk::ImageFileWriter<rawImageType_uint> labelWriterType;
typedef itk::ImageFileWriter<rawImageType_8bit> rawWriterType_8;

typedef itk::RescaleIntensityImageFilter< rawImageType_16bit, rawImageType_8bit > RescaleFilterType;
typedef itk::RescaleIntensityImageFilter< rawImageType_16bit, rawImageType_16bit > RescaleFilterType_16to16;
typedef itk::RescaleIntensityImageFilter< rawImageType_flo, rawImageType_flo > RescaleFilterType_floToflo;

typedef itk::CastImageFilter< rawImageType_16bit, rawImageType_flo > CasterFilterType_16Tofloat;

typedef itk::MedianImageFilter<rawImageType_flo, rawImageType_flo> MedianFilterType_floToflo;

typedef itk::ImageRegionIterator< rawImageType_8bit > IteratorType_8bit;
typedef itk::ImageRegionIterator< rawImageType_16bit > IteratorType_16bit;
typedef itk::ImageRegionIterator< rawImageType_uint > IteratorType_uint;


typedef itk::RGBPixel<unsigned char> RGBPixelType_8bit;
typedef itk::Image<RGBPixelType_8bit,3> RGBImageType_8bit;
typedef itk::ImageRegionIterator<RGBImageType_8bit> RGBIteratorType_8bit;