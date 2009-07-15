// File: tracing_config.h
//Code generated from 55 matlab m files
#ifndef TRACECONFIG_H
#define TRACECONFIG_H

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "itkImage.h"

//#include <libxml/parser.h>
//#include <libxml/tree.h>
//#include <libxml/xmlreader.h>


#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

typedef itk::Image <float, 2> ImageType2D;
typedef itk::Image <float, 3> ImageType3D;
typedef vnl_vector_fixed<double,3> Vect3;
typedef vnl_matrix_fixed <double,3,3> Mat33;


#define DEBUG_LEVEL 1



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
