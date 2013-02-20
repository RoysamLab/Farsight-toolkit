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
 \brief This file contains the main function of the testing class. The user can input a 3D image and a seed point [x,y,z] coordinate location. The tracer will start from the specifiled location and travel in both directions to output a trace. 
 \note The image should have intensities between 0 and 255, with the foreground dark (closer to 0) and background white (closer to 1)
 \author $ Author: Amit Mukherjee, James Alex Tyrrell $
 \version $Revision 1.0 $
 \date May 05 2009
*/

#include <iostream>
#include <fstream>

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkStatisticsImageFilter.h"
#include <vnl/vnl_vector_fixed.h>

#include "TraceContainer3D.h"
#include "NodeContainer3D.h"
#include "TVessel.h"
#include "SegInit.h"
#include "Trace.h"


#define STEP_FWD 1
#define STEP_BWD -1

typedef itk::Image< float, 2 >   ImageType2D;
typedef itk::Image< float, 3 >   ImageType3D;

typedef vnl_vector_fixed<double,3> Vect3;



int main (int argc, char *argv[])
{

	if (argc != 5)	{
		std::cout << "Usage: "<<argv[0] << " [ImageFile] [Seedx] [Seedy] [Seedz]" <<std::endl;
		return EXIT_FAILURE;
	}

	std::string fname(argv[1]);

	std::cout << "Tracing "<< fname.c_str() << std::endl;

	ImageType3D::Pointer image3D = ImageType3D::New();
	ImageType2D::Pointer MIPimage = ImageType2D::New();

	typedef itk::ImageFileReader<ImageType3D> ReaderType;
	ReaderType::GlobalWarningDisplayOff();
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(fname.c_str() );
	image3D = (reader->GetOutput());
	image3D->Update();
	std::cout << "Image of size " << image3D->GetBufferedRegion().GetSize() << " read successfully " << std::endl;
	
	try	{
		long segID = 0;

		TVessel* seg = new  TVessel();
		seg->mu[0] = atof(argv[2]);
		seg->mu[1] = atof(argv[3]);
		seg->mu[2] = atof(argv[4]);
		seg->ID = segID;	

		SegInit *fitter = new SegInit();
		bool ret = fitter->fitSE(image3D, *seg, 100, 1.35);
		delete fitter;

		NodeContainer3D::Pointer NodeList = NodeContainer3D::New();
		NodeList->AddNode(seg);

		Trace * tr = new Trace();
		tr->TraceID = 1;	// tID set as 1 arbitrarily
		tr->numNodes = 1;
		tr->L = seg->L;
		tr->NodeAID = seg->ID;
		tr->NodeBID = seg->ID;

		for (int i=0; i<3; i++)	{
			tr->dirA[i] = seg->R3[i];
			tr->dirB[i] = -1*seg->R3[i];
		}
		tr->PrintSelf();

		seg->TraceID = 1;

		seg->PrintSelf();
		int numStep = 0;
		segID = 1;
		std::cout <<"Starting Iteration"<<std::endl;

		std::string SegTxtFname = fname + std::string("Output.txt");

		while(numStep < 1000)	{

			numStep++;
			TVessel *seg1 = tr->Step(seg, image3D, segID++, STEP_FWD);
			if (!seg1) {
				std::cout <<"seg 1 returned empty"<< std::endl;
				break;
			}

			if (seg1->IsSegmentValid(seg)==0)	{
				seg1->PrintSelf();
				std::cout <<"Seg Valid Violation"<< std::endl;
				break;
			}


			NodeList->AddNode(seg1);

			tr->UpdateTrace(seg, seg1, 1);
			NodeList->getPropagationDirectionR3(seg1->ID, tr);


			seg1->PrintSelf();	tr->PrintSelf();	NodeList->PrintSelf();	std::cout <<"Step done "<< std::endl;


			if (NodeList->IsTraceValid(tr, seg1)==0)	{
				std::cout <<"Trace Valid Violation"<< std::endl;
				break;
			}

			long hitSeg;
			if (hitSeg = NodeList->HitTest(seg1, tr) > -1)	{
				std::cout <<"Hit segment" << hitSeg<< std::endl;
				break;
			}

			seg = seg1;
		}


		//travel reverse in node container
		std::cout <<"Starting reverse tracing"<< std::endl;
		seg = NodeList->getSegment(0);
		tr->Reverse( seg );

		numStep = 0;
		while(numStep < 1000)	{
			numStep++;
			TVessel *seg1 = tr->Step(seg, image3D, segID++, STEP_BWD);

			if (!seg1) {
				std::cout <<"seg 1 returned empty"<< std::endl;
				break;
			}

			if (seg1->IsSegmentValid(seg)==0)	{
				seg1->PrintSelf();
				std::cout <<"Seg Valid Violation"<< std::endl;
				break;
			}

			NodeList->AddFront(seg1);
			tr->UpdateTrace(seg, seg1, 2);
			NodeList->getPropagationDirectionR3(seg1->ID, tr);


			seg1->PrintSelf();	tr->PrintSelf();	NodeList->PrintSelf();	std::cout <<"Step done "<< std::endl;


			if (NodeList->IsTraceValid(tr, seg1)==0)	{
				std::cout <<"Trace Valid Violation"<< std::endl;
				break;
			}

			long hitSeg;
			if (hitSeg = NodeList->HitTest(seg1, tr) > -1)	{
				std::cout <<"Hit segment# " << hitSeg<< std::endl;
				break;
			}

			seg = seg1;
		}

		NodeList->WriteSegmentsToTextFile(SegTxtFname);
		NodeList->WriteSegmentsToXMLFile(SegTxtFname);

		delete tr;

	}
	catch (itk::ExceptionObject & err )	    {
	    std::cout << "ExceptionObject caught !" << std::endl;
	    std::cout << err << std::endl;
	    return EXIT_FAILURE;
    }

	return EXIT_SUCCESS;
}


