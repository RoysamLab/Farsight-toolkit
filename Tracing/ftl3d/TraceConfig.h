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

/**
 \brief The Stores configuration info of the tracing. 
 \author $ Author:  Amit Mukherjee $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit , Rensselaer Polytechnic institute Troy NY 12180.

#ifndef TRACECONFIG_H
#define TRACECONFIG_H

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "itkImage.h"
#include "itkFixedArray.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "tinyxml.h"

typedef itk::Image <float, 2> ImageType2D;
typedef itk::Image <float, 3> ImageType3D;
typedef vnl_vector_fixed<double,3> Vect3;
typedef vnl_matrix_fixed <double,3,3> Mat33;


class TraceConfig: public itk::LightObject        {

public:
	   typedef TraceConfig Self;
	   typedef  itk::SmartPointer< Self >  Pointer;
	   typedef  itk::SmartPointer< const Self >  ConstPointer;
	   itkNewMacro(Self);
	   //default constructor loads default configuration values

		bool LoadParameters (char *);

		int getGridSpacing() {return(GridSpacing);}
		double getStepRatio() {return(StepRatio);}
		double getAspectRatio() {return(AspectRatio);}
		double getTHRESHOLD() {return THRESHOLD;}
		double getminContrast() {return minContrast;}
		double getMaximumVesselWidth() {return(MaximumVesselWidth);}
		double getMinimumVesselLength() {return(MinimumVesselLength);}
		int getFitIterations() { return FitIterations;}
		double getMinimumVesselWidth() {return(MinimumVesselWidth);}
		double getStartThreshold() {return StartTHRESHOLD;}
		itk::FixedArray<double,3> GetSpacing() {return Spacing;}
		int getNumberOfDataFiles() {return numDataFiles;}
		int getHessianFlag() {return UseMultiscaleHessianFilter;}
		
		std::string getInputFileName(unsigned int i) {
			if (i < numDataFiles)
				return(InputFileNames[i]);
			else
				return (std::string(" "));
		}
		std::string getOutputFileName(unsigned int i) {
			if (i < numDataFiles)
				return(OutputFileNames[i]);
			else
				return (std::string(" "));
		}
		
		void SetFileNames(char* fname);
		void SetGridSpacing(char *);
		void SetAspectRatio(char *);

		~TraceConfig ();

private:
		bool ParseDoubleInput(double, double, double, const char* );
		std::vector<std::string> InputFileNames;
        std::vector<std::string> OutputFileNames;
		unsigned int numDataFiles;

        //tracing parameters
		int GridSpacing;
		double StepRatio;
		double AspectRatio;
		double THRESHOLD;
		double minContrast;
		double MaximumVesselWidth;
		double MinimumVesselLength;
		int FitIterations;
		double MinimumVesselWidth;
		double StartTHRESHOLD;
		itk::FixedArray<double, 3> Spacing;
		unsigned int UseMultiscaleHessianFilter;

protected:
   		TraceConfig ();

};


#endif
