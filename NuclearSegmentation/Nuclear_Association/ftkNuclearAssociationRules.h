/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Module:    $RCSfile: ftkNuclearAssociationRules.h,v $
  Language:  C++
  Date:      $Date: 2008/11/27 1:00:00 $
  Version:   $Revision: 1 $
 
=========================================================================*/
#ifndef __ftkNuclearAssociationRules_h
#define __ftkNuclearAssociationRules_h

#include <limits.h>

#include <ftkObject.h>
#include <ftkFeatures/ftkObjectAssociation.h>
#include "ftkImage.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageDuplicator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itkSignedDanielssonDistanceMapImageFilter.h>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <vtkTable.h>

#include "VolumeOfInterest.h"

#ifdef _OPENMP
#include "omp.h"
#endif

typedef itk::Image< unsigned short, 3 > LabImageType;
typedef itk::Image< unsigned short, 3 > TargImageType;
std::vector<float> compute_ec_features( TargImageType::Pointer input_image,  LabImageType::Pointer input_labeled, int number_of_rois, unsigned short thresh, int surr_dist, int inside_dist  );
unsigned short returnthresh( TargImageType::Pointer input_image, int num_bin_levs, int num_in_fg );

namespace ftk
{ 

class AssociativeFeatureCalculator
{
public:
	AssociativeFeatureCalculator();
	~AssociativeFeatureCalculator();
	void SetInputFile(std::string filename);
	void SetInputs(ftk::Image::Pointer inp_labeled_image, int inp_channel_number, ftk::Image::Pointer seg_labeled_image, int seg_channel_number, ftk::AssociationRule *associationrule);
	void SetFeaturePrefix(std::string prefix);			//Set Prefix for feature names
	vtkSmartPointer<vtkTable> Compute();			//Compute and return table with values (for all objects)
	void Update(vtkSmartPointer<vtkTable> table);		//Update the features in this table whose names match (sets doFeat)
	void Append(vtkSmartPointer<vtkTable> table);		//Compute features that are ON and append them to the existing table
	bool AreInputsSet(){ return inputs_set; }

private:
	ftk::AssociationRule *input_association;
	std::string inFilename;
	std::string fPrefix;
	bool inputs_set;
	LabImageType::Pointer lab_im;
	TargImageType::Pointer inp_im;
};


/** \class NuclearAssociation
 *  \brief To define a set of association rules between nuclei and other objects and to compute the corresponding associative measures. 
 *  
 * \author Yousef Al-Kofahi, Rensselear Polytechnic Institute (RPI), Troy, NY USA
 */
class NuclearAssociationRules : public ObjectAssociation
{
public:
	/* Contsructor */
	NuclearAssociationRules(std::string segImageName, int numOfAssocRules);
	NuclearAssociationRules(std::string AssocFName, int numOfRules, LabImageType::Pointer lImage, TargImageType::Pointer iImage);
	~NuclearAssociationRules();

	/* This method computes all the associative measurements for all the objects */
	void Compute(vtkSmartPointer<vtkTable> tbl = NULL);
		
	/* Get the number of objects */
	int GetNumOfObjects() {return numOfLabels;};
private:
	/* Private member variables */
	typedef itk::Image< double, 3 > DistImageType;
	typedef itk::ImageFileReader< LabImageType > ReaderType;
	typedef itk::ImageFileWriter< LabImageType > WriterType;
	typedef itk::BinaryThresholdImageFilter< LabImageType, LabImageType > BinaryThresholdType;
	typedef itk::LabelGeometryImageFilter< LabImageType, LabImageType > LabelGeometryType;
	//typedef itk::ApproximateSignedDistanceMapImageFilter<DistImageType, DistImageType > DTFilter ;
	typedef itk::SignedDanielssonDistanceMapImageFilter<DistImageType, DistImageType > DTFilter ;

	LabImageType::Pointer labImage;
	TargImageType::Pointer inpImage;
	LabelGeometryType::Pointer labGeometryFilter;	
	int x_Size;
	int y_Size;
	int z_Size;
	int imDim;
	bool inputs_set;

	unsigned short thresh;


private:
	/* This method is used to compute a single associative measurement for one cell */
	float ComputeOneAssocMeasurement(itk::SmartPointer<TargImageType> trgIm, int ruleID, int objID);

	/* the next functions will compute the measurements for the different cases */
	float FindMin(std::vector<int> LST);
	float FindMax(std::vector<int> LST);
	float ComputeTotal(std::vector<int> LST);
	float ComputeAverage(std::vector<int> LST);

	int num_rois;
	
}; // end NuclearAssociation
}  // end namespace ftk

#endif	// end __ftkNuclearAssociationRules_h
