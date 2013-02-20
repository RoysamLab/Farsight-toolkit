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
#include "mdlVolumeProcess.h"
#include "mdlIntegratedSkeleton.h"
#include "mdlMST.h"
#include "mdlBSplineFitting.h"
#include "mdlUtils.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char *argv[])
{
	if(argc < 3)
    {
		std::cerr << "Usage: " << argv[0] << " inputImageFile  inputVtkFile  outputVtkFile " <<std::endl;
		return EXIT_FAILURE;
    }
    std::string inputImageFile = argv[1];
	std::string inputVtkFile = argv[2];
	std::string outputVtkFile = argv[3];
	

	// IMAGE READER
	typedef itk::ImageFileReader< mdl::ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(inputImageFile);
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "READER FAILED: " << err << std::endl ;
		return EXIT_FAILURE;
	}
	mdl::ImageType::Pointer img = reader->GetOutput();
	
 
	
	mdl::vtkFileHandler * fhd0 = new mdl::vtkFileHandler();
    fhd0->GetNodesandLinesFromVtkfile(inputVtkFile);
    std::vector<mdl::fPoint3D> nodes = fhd0->getNodes();
    std::vector<mdl::pairE> bbpairs =  fhd0->getLines();

	std::cerr << "BSPLINE\n";
	
	mdl::BSplineFitting *bspline = new mdl::BSplineFitting( img );
	bspline->SetDebug(false);
	bspline->SetLevels(6);
	bspline->SetOrder(3);
	bspline->SetNodes( &nodes );
	bspline->SetBBPairs( &bbpairs );
	bspline->Update();
	nodes = bspline->GetNodes();
	bbpairs = bspline->GetBBPairs();
	delete bspline;


	fhd0->SetNodes(&nodes);
	fhd0->SetLines(&bbpairs);
	fhd0->Write(outputVtkFile.c_str());
	delete fhd0;

	return EXIT_SUCCESS;

}
