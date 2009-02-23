/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef __ftkFeatures_h
#define __ftkFeatures_h

#include <map>

namespace ftk
{

//********************************
// USE UNION TO HOLD THE FEATURE DATA
//********************************
typedef union {
	float f;
	int i;
} Value;

typedef struct {
	bool isInteger;
	std::string units;
	std::string description;
} Info;

//*****************************************************************************************
// THIS MAP HOLDS ALL FEATURE INFORMATION FOR ONE OBJECT
//*****************************************************************************************
typedef std::map< std::string, Value > LabelImageFeatureValueMapType;
typedef std::map< std::string, Info > LabelImageFeatureInfoMapType;
 

}  // end namespace ftk

#endif