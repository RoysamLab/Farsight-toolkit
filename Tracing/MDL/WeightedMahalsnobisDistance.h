/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include <stdio.h>
#include <cstring>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "itkArray.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkWeightedMeanCalculator.h"
#include "itkWeightedCovarianceCalculator.h"
// Software Guide : EndCodeSnippet
using std::cerr;
using std::cout;
using std::cin;
using std::endl;

// #define FeatureNumber  3

typedef itk::Vector< float, 3 > MeasurementVectorType;

class ExampleWeightFunction :
  public itk::FunctionBase< MeasurementVectorType, double >
{
public:
  /** Standard class typedefs. */
  typedef ExampleWeightFunction Self;
  typedef itk::FunctionBase< MeasurementVectorType, double > Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  
  /** Standard macros. */
  itkTypeMacro(ExampleWeightFunction, FunctionBase);
  itkNewMacro(Self);

  /** Input type */
  typedef MeasurementVectorType InputType;

  /** Output type */
  typedef double OutputType;

  /**Evaluate at the specified input position */
  OutputType Evaluate( const InputType& input ) const 
    {
      if ( input[0] < 3.0 )
        {
        return 0.5;
        }
      else
        {
        return 0.01;
        }
    }

protected:
  ExampleWeightFunction() {}
  ~ExampleWeightFunction() {}
} ; // end of class


class MDLClassifier
{ 
   int FeatureNumber;
   double *MeanVector; 
   double *CovarianceMatrix; 
   double *InverseCovarianceMatrix;
  
public:
   
  MDLClassifier(int a)
   {
		FeatureNumber = a;
		if (FeatureNumber>0)
		{
	    //Sample = new double [a];
		MeanVector = new double [FeatureNumber];
        CovarianceMatrix =  new double [FeatureNumber*FeatureNumber];
        InverseCovarianceMatrix = new double [FeatureNumber*FeatureNumber];
		}
   }
    int MeanVectorandVarianceMatrix(char *filename);
    double MahalanobisDist(double *sample); 
    double MahalanobisDist(double meanDensityBranch, double length_leaf, double meanVesselBranch, int spineOne);
    void print()
    {
     int i,j;
     for (i=0;i<FeatureNumber;i++)
	 std::cout << " "  << MeanVector[i];
     std::cout << std::endl;
     int IDX =0;

     for (i=0;i<FeatureNumber;i++)
     for (j=0;j<FeatureNumber;j++)
     {   IDX= i*FeatureNumber + j;
		 std::cout << " "  << *(InverseCovarianceMatrix+IDX);
     }

    }

  ~MDLClassifier(void)
  {
	  if (FeatureNumber >0)
	  {
		  //delete []Sample;
	      delete []MeanVector;
		  delete []CovarianceMatrix;
		  delete []InverseCovarianceMatrix ;
	  }
  }
}; // end of class 