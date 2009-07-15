/////////////////////////////////////////////////////////////////////////////////
// FILE: TraceMain.cpp
//

#include <iostream>
#include <fstream>

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkStatisticsImageFilter.h"
#include "vnl\vnl_vector_fixed.h"

#include "TraceContainer3D.h"
#include "NodeContainer3D.h"
#include "TVessel.h"
#include "SegInit.h"
#include "Trace.h"

#include "myDebug.h"

#define STEP_FWD 1
#define STEP_BWD -1

typedef itk::Image< float, 2 >   ImageType2D;
typedef itk::Image< float, 3 >   ImageType3D;

typedef vnl_vector_fixed<double,3> Vect3;

void GenerateMIP(ImageType3D::Pointer, std::string&);


int main (int argc, char *argv[])
{

	if (argc != 2)	{
		std::cout << "Usage: "<<argv[0] << " [ImageFile]" <<std::endl;
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
	std::string MIPFname("MIPTestImage.png");
    //GenerateMIP( image3D, MIPFname );

	try	{
		long segID = 0;

		TVessel* seg = new  TVessel();
		seg->mu[0] = 20.0;
		seg->mu[1] = 11.0;
		seg->mu[2] = 21.0;
		seg->f = 120;
		seg->ID = segID;	// ID set as 1

		SegInit *fitter = new SegInit();
		bool ret = fitter->fitSE(image3D, *seg, 100, 1.35);
		delete fitter;

		NodeContainer3D::Pointer NodeList = NodeContainer3D::New();
		NodeList->AddNode(seg);

		Trace * tr = new Trace();
		tr->TraceID = 1;	// tID set as 1
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
		S("Starting Iteration")

		std::string SegTxtFname("SegAllFile.txt");

		while(numStep < 1000)	{

			numStep++;
			TVessel *seg1 = tr->Step(seg, image3D, segID++, STEP_FWD);
			if (!seg1) {
				S("seg 1 returned empty")
				break;
			}

			if (seg1->IsSegmentValid(seg)==0)	{
				seg1->PrintSelf();
				S("Seg Valid Violation")
				break;
			}


			NodeList->AddNode(seg1);

			tr->UpdateTrace(seg, seg1, 1);
			NodeList->getPropagationDirectionR3(seg1->ID, tr);


			seg1->PrintSelf();	tr->PrintSelf();	NodeList->PrintSelf();	S("Step done ")


			if (NodeList->IsTraceValid(tr, seg1)==0)	{
				S("Trace Valid Violation")
				break;
			}

			long hitSeg;
			if (hitSeg = NodeList->HitTest(seg1, tr) > -1)	{
				S("Hit segment" << hitSeg)
				break;
			}

			seg = seg1;
		}


		//travel reverse in node container
		S("Starting reverse tracing")
		seg = NodeList->getSegment(0);
		tr->Reverse( seg );

		numStep = 0;
		while(numStep < 1000)	{
			numStep++;
			TVessel *seg1 = tr->Step(seg, image3D, segID++, STEP_BWD);

			if (!seg1) {
				S("seg 1 returned empty")
				break;
			}

			if (seg1->IsSegmentValid(seg)==0)	{
				seg1->PrintSelf();
				S("Seg Valid Violation")
				break;
			}


			NodeList->AddFront(seg1);

			tr->UpdateTrace(seg, seg1, 2);
			NodeList->getPropagationDirectionR3(seg1->ID, tr);


			seg1->PrintSelf();	tr->PrintSelf();	NodeList->PrintSelf();	S("Step done ")


			if (NodeList->IsTraceValid(tr, seg1)==0)	{
				S("Trace Valid Violation")
				break;
			}

			long hitSeg;
			if (hitSeg = NodeList->HitTest(seg1, tr) > -1)	{
				S("Hit segment" << hitSeg)
				//break;
			}

			seg = seg1;
		}

		NodeList->WriteSegmentsToTextFile(SegTxtFname);
		delete tr;

	}
	catch (itk::ExceptionObject & err )	    {
	    std::cout << "ExceptionObject caught !" << std::endl;
	    std::cout << err << std::endl;
	    return EXIT_FAILURE;
    }

	return EXIT_SUCCESS;
}



void GenerateMIP(ImageType3D::Pointer im3D, std::string& MIPfile)	{
	ImageType3D::SizeType sz3 = im3D->GetRequestedRegion().GetSize();
	typedef itk::Image<unsigned char,2> CharImageType2D;

	CharImageType2D::SizeType sz2 = {{sz3[0], sz3[1] }};

	P(sz3)
	CharImageType2D::Pointer im2D = CharImageType2D::New();
	im2D->SetRegions(sz2 );
	im2D->Allocate();
	P(sz3)

	for (unsigned long x=0; x<sz3[0]; x++)	{
		for (unsigned long y=0; y<sz3[1]; y++)	{

			//ImageType3D::PixelType maxVal = 0.0;
			ImageType3D::PixelType minVal = 255.0;

			for (unsigned long z=0; z<sz3[2]; z++)	{
				ImageType3D::IndexType nd3 = {{x,y,z}};
				ImageType3D::PixelType val = im3D->GetPixel(nd3);
			//	maxVal = (maxVal > val) ? maxVal : val;
				minVal = (minVal < val) ? minVal : val;
			}

			CharImageType2D::IndexType nd2 = {{x, y}};
			//im2D->SetPixel(nd2, maxVal);
			im2D->SetPixel(nd2, static_cast<unsigned char> (minVal));
		}
	}

	itk::ImageFileWriter<CharImageType2D>::Pointer wt = itk::ImageFileWriter<CharImageType2D>::New();
	wt->SetFileName(MIPfile.c_str());
	wt->SetInput(im2D);
	wt->Update();
}
