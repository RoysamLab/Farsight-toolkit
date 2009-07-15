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

// File: tracing_config.h
//Code generated from 55 matlab m files
#ifndef TRACECONFIG_H
#define TRACECONFIG_H

#include <itkImage.h>
#include <itkSmartPointer.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlreader.h>


#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

typedef itk::Image< float, 2 >   ImageType;
typedef vnl_vector_fixed <double,2> Vect2;
typedef vnl_matrix_fixed <double,2,2> Mat22;


#define DEBUG_LEVEL 1

/*
DEBUG_LEVEL		0  	- No message
DEBUG_LEVEL		1	- Few
DEBUG_LEVEL		2	- Detials
*/

class TraceConfig        {
public:

	   //default constructor loads default configuration values
	   TraceConfig ();
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
	   double getStepSize(){return(StepSize);}
	  /* void setGuiParameters(double GS, double SR, double AR, double Sensitivity)	{
		   this->GridSpacing = (int)GS;
		   this->StepRatio = SR;
		   this->AspectRatio = AS;
		   this->PROP = Sensitivity;
	   }*/
		
	   std::string getInputFileName() {return(InputFileName);}
	   std::string getOutputFileName() {return(OutputFileName);}



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
        double StepSize;

};


#endif
