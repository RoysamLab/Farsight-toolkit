#ifndef _MICROGLIA_MOVIE_SEG_TRACER_H
#define _MICROGLIA_MOVIE_SEG_TRACER_H

#include <itkImage.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "..\Tracing\MultipleNeuronTracer\MultipleNeuronTracer.h"

class MicrogliaMovieSegTracer
{
public: 
	MicrogliaMovieSegTracer();
	virtual ~MicrogliaMovieSegTracer();
	void StartSegTracing(const char* filename, const char* paramFileName);

	typedef MultipleNeuronTracer::CharImageType3D CharImageType;
	typedef MultipleNeuronTracer::ImageType3D FloatImageType;
	typedef MultipleNeuronTracer::LabelImageType3D LabelImageType;
	typedef itk::CastImageFilter< CharImageType, FloatImageType> CasterTypeCharToFloat;
	typedef itk::ImageFileReader< CharImageType> ReaderType;
	typedef itk::ImageFileWriter< CharImageType> CharWriterType;
	typedef itk::CastImageFilter< LabelImageType, CharImageType> CasterTypeLabelToChar;

protected:
	void ReadXMLProjectFile(const char* filename);
	void LoadOptions(const char* paramFileName);
	CharImageType::Pointer ReadImage(const char * fileName);
	void WriteImage(const char* filename, LabelImageType::Pointer labelImage);
	template <class T> bool SetParamValue(std::map<std::string,std::string> &opts, std::string str, T &value, T defVal);

	void GenerateSeedPointsForSoma(CharImageType::Pointer inputImage, std::vector< itk::Index<3> > &seedVec);
	LabelImageType::Pointer SegmentSoma(FloatImageType::Pointer inputImage, std::vector< itk::Index<3> > &seedVec);
	void StartTracing(FloatImageType::Pointer inputImage, std::vector< itk::Index<3> > &seedVec, LabelImageType::Pointer somaImage);
	void CalculateLMeasures();

private:
	std::vector< std::string> files;

	/// For SeedDetection
	int shift;
	double scaleMin;
	double scaleMax;
	double regionXY;
	double regionZ;
	int useDistMap;
	int sampling_ratio_XY_to_Z;

	/// For Soma
	double alfa;
	double beta;
	int timethreshold; 
	double curvatureScaling; 
	double rmsThres;
	int holeSize;
	int minObjSize;

	/// For Tracing
	float intensity_threshold;
	float contrast_threshold;
	int cost_threshold;
	float debris_threshold;
	int offshoot;
	int device;
};

#endif
