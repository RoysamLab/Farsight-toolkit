#ifndef CURVELET_H_
#define CURVELET_H_
/*
Copyright (C) 2004 Caltech
Written by Lexing Ying
*/

#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"

#include <stdlib.h>
#include <time.h>
#include <omp.h>

//itk includes
#include "itkImage.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>

//typedefs
typedef unsigned char InputPixelType;
typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<InputPixelType,2> Input2DImageType;
typedef itk::Image<float,3> FloatImageType;
typedef itk::Image<float,2> Float2DImageType;


#define MIN(a,b) (((a) > (b))? (b) : (a))
#define MAX(a,b) (((a) < (b))? (b) : (a))

class Curvelet
{
public:
	Curvelet();
	~Curvelet();
	Curvelet(InputImageType::Pointer NewInputImage);
	InputImageType::Pointer RunOnInputImage(InputImageType::Pointer InputImage);
	void SetSigma(float sigma);

	template<typename T>
	void circshift(T input, int shiftx, int shifty, T &output);
	Input2DImageType::Pointer getSlice(InputImageType::Pointer im, int slice);
	template <typename PixelType>
	void copyslice(typename itk::Image<PixelType,2>::Pointer im1, typename itk::Image<PixelType,3>::Pointer im2, int slice);
	void getCurveletsForOneSlice(Input2DImageType::Pointer im,Input2DImageType::Pointer &om,Float2DImageType::Pointer &cosim, Float2DImageType::Pointer &sinim);


private:
	int nbangles_coarse;
	int ac;
	int nshifts;
	float neighb_weight;
	float nsigmas_coarse;
	float nsigmas_fine;
	float tuning_neighb;
	float sigma_ratio;
	int numt;
	int tile_size;
	int border;
	int slices;
	//InputImageType::Pointer InputImage;
};
#endif