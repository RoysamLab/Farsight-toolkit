#ifndef __GradientEnhancedDistanceMapImageFilter_H
#define __GradientEnhancedDistanceMapImageFilter_H



#include <iostream>
#include <algorithm>
#include <float.h>
#include <math.h>
#include <queue>

//ITK include
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGradientMagnitudeImageFilter.h"


struct im_ind{
	int i;
	int j;
};

int gradient_enhanced_distance_map_2D( float *GRAD_IM_2D, float *GRAD_IMW, int size1, int size2 );
//added by Yousef on 09/14/2009
int computeGradientImage(unsigned char *IM, float* gradIM, int r, int c, int z);
void prepareInputImage(unsigned short *binImage, unsigned short* seedsImage, float* outImage, int r, int c, int z);

#endif