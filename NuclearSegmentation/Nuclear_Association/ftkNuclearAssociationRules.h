/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Module:    $RCSfile: ftkNuclearAssociationRules.h,v $
  Language:  C++
  Date:      $Date: 2008/11/27 1:00:00 $
  Version:   $Revision: 1 $
 
=========================================================================*/
#ifndef __ftkNuclearAssociationRules_h
#define __ftkNuclearAssociationRules_h

#include <ftkObject.h>
#include <SegmentationCommon/ftkObjectAssociation.h>
#include "itkLabelGeometryImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkApproximateSignedDistanceMapImageFilter.h"

namespace ftk
{ 
/** \class NuclearAssociation
 *  \brief To define a set of association rules between nuclei and other objects and to compute the corresponding associative measures. 
 *  
 * \author Yousef Al-Kofahi, Rensselear Polytechnic Institute (RPI), Troy, NY USA
 */
class NuclearAssociationRules : public ObjectAssociation
{
public:
	/* Contsructor */
	NuclearAssociationRules(string segImageName, int numOfAssocRules);	

	/* This method computes all the associative measurements for all the objects */
	void Compute();
		
	/* Get the number of objects */
	int GetNumOfObjects() {return numOfLabels;};
private:
	/* Private member variables */
	typedef itk::Image< unsigned short, 3 > LabImageType;
	typedef itk::Image< unsigned short, 3 > TargImageType;
	typedef itk::Image< float, 3 > DistImageType;
	typedef itk::ImageFileReader< LabImageType > ReaderType;
	typedef itk::LabelGeometryImageFilter< LabImageType, LabImageType > LabelGeometryType;
	typedef itk::ApproximateSignedDistanceMapImageFilter<DistImageType, DistImageType > DTFilter ;
	LabImageType::Pointer labImage;
	LabelGeometryType::Pointer labGeometryFilter;	
	int x_Size;
	int y_Size;
	int z_Size;
	int imDim;		


private:
	/* This method is used to compute a single associative measurement for one cell */
	float ComputeOneAssocMeasurement(itk::SmartPointer<TargImageType> trgIm, int ruleID, int objID);

	/* the next functions will compute the measurements for the different cases */
	float FindMin(vector<int> LST);
	float FindMax(vector<int> LST);
	float ComputeTotal(vector<int> LST);
	float ComputeAverage(vector<int> LST);
	
}; // end NuclearAssociation
}  // end namespace ftk

#endif	// end __ftkNuclearAssociationRules_h
