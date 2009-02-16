/** @file nucleus_feature_util.h
*   @brief containing utility functions that compute features for a nucleus
*
*   @author Chia-Ling Tsai (2007/09/28)
*/

#ifndef _NUCLEUS_FEATURE_UTIL_H_
#define _NUCLEUS_FEATURE_UTIL_H_

#include "itkImage.h"
#include "rich_cell.h"

typedef unsigned char                        IntensityPixelType;
typedef itk::Image< IntensityPixelType, 3 >  IntensityImageType;
typedef float                                InternalPixelType;
typedef itk::Image< InternalPixelType, 3 >   InternalImageType;
typedef itk::Image< unsigned short,3 >       LabelImageType;

//: Compute features that are related to the volume only
//
//  Features comptued are: 
void
compute_volume_related_features( rich_cell::Pointer cell );

//: Compute features related to intensity information of a cell
//
//  Features computed are: 1 ave_intensity_ 2 avg_bound_intensity_, 3 
void 
compute_intensity_related_features( IntensityImageType::Pointer image, rich_cell::Pointer cell);

//: Compute features related to the neighbors
//
//  Features computed are:
void 
compute_neighbor_related_features( LabelImageType::Pointer mask, rich_cell::Pointer cell);

#endif
