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

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

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

	   void LoadParameters (char *);
	   int getGridSpacing() {return(GridSpacing);}
	   int getInitSeedSize() {return(SeedSize);}
	   int getSeedIntensityThreshold() {return(SeedIntensityThreshold);}
	   double getStepRatio() {return(StepRatio);}
	   double getAspectRatio() {return(AspectRatio);}
	   double getMinimumVesselWidth() {return(MinimumVesselWidth);}
	   double getMaximumVesselWidth() {return(MaximumVesselWidth);}
	   double getMinimumVesselLength() {return(MinimumVesselLength);}
	   double getPROP() {return(PROP);}

	   std::string getInputFileName() {return(InputFileName);}
	   std::string getOutputFileName() {return(OutputFileName);}
	   void SetFileNames(char* fname);
	   void SetGridSpacing(char *);
	   void SetAspectRatio(char *);

		~TraceConfig ();

private:
        std::string InputFileName;
        std::string OutputFileName;

		//seed parameters
        int GridSpacing;
        int SeedSize;
        int SeedIntensityThreshold;
        //tracing parameters
        double StepRatio;
        double AspectRatio;
        double MinimumVesselWidth;
        double MaximumVesselWidth;
        double MinimumVesselLength;
        double PROP;


protected:
   		TraceConfig ();

};


#endif
