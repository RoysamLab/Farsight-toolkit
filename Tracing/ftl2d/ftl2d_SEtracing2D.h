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

/** @file SEtracing2D.cpp
*   @brief The main program of the tracing algorithm using Superellipsoid.
*
*   @author Amit Mukherjee
*////////////////////////////////////////////////////////////////////////////////// 

#ifndef SETRACING2D_H
#define SETRACING2D_H

#include "Common/farsight_object.h"

#include "ftl2d_TraceConfig.h"
#include "ftl2d_SeedContainer2D.h"
#include "ftl2d_Tracer.h"

#include <itkImage.h>
#include <itkSmartPointer.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

typedef itk::Image<unsigned char, 2>::ConstPointer LabelImageTypePointer;


class SETracing2D: public farsight_object<2>	{
	public:
		SETracing2D() {}
	  
		//: destructor
  		virtual ~SETracing2D(){}
  		
  		typedef SETracing2D			Self;
  		typedef farsight_object<2> SuperClass;
  		
  		
  		typedef itk::Image <float, 2> ImageType;
  		typedef itk::Image<unsigned char,2> LabelImageType;
		typedef itk::ImageFileReader< ImageType >  ReaderType;
	//	typedef SuperClass::ImageConstPtrType LabelImageTypePointer;
  		
  		//reads the tracing parameter 
  		virtual void read_xml(std::string const& file_path, std::string const& xml_filename);

		// set the image name and path as well as loads the image to be traced
  		void set_image( std::string const& image_path, std::string const& image_name);
  		
  		// main tracing module
  		bool run();
  		
  		//writes the tracing output 
  		void write_xml(std::string const& file_path, std::string const& xml_filename);
  		
  		//creates a unsigned char image and return 1 in places where traces exist and 0 at background 
  		// for overlaying on the trace. 
  		virtual ImageConstPtrType display() ;
  		
  		void WriteTraceLabel();
  		
	
	
	private:
		TraceConfig * m_Config;					// config
		SeedContainer2D * m_Seed;				// seed
		Tracer * m_Tracer;						// tracer
		ImageType::Pointer m_image;  			// image data
		LabelImageType::Pointer m_Label;		// label image
		std::string		m_LabelImageName;
		
	/*protected:
	  	std::vector< std::string > parameters_;
	  	ImagePtrType               label_image_;
	  	std::string                image_name_;
  	  	std::string                image_path_;
  */
};
	
#endif
