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
#include "ftkPreprocess2.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <string.h>
#include <tinyxml/tinyxml.h>

void usage(const char *funcName);
void pipeline(std::string pipeName, ftk::Preprocess * prep);

int main(int argc, char *argv[])
{
	if( argc < 3 )
	{
		usage(argv[0]);
		std::cerr << "PRESS ENTER TO EXIT\n";
		getchar();
		return EXIT_FAILURE;
	}
	
	std::string inputFilename = "";					// Name of the input image;
	std::string filterName = "";					// Name of the filter to apply
	std::vector<double> params;						// Params for the filter to apply
	std::string pipeName = "";						// Filename of pipeline file
	std::string outputFilename = "prep_out.tif";	// Name of the output image;
	bool useColor = false;
	char color = 'w';								// Create intensity image from colors

	int c = 1;
	std::string token = argv[1];
	bool done = false;
	while(!done)
	{
		if( token == "-i" )
		{
			if(c+1 < argc)
				inputFilename = argv[++c];
		}
		else if( token == "-p" )
		{
			if(c+1 < argc)
				pipeName = argv[++c];
		}
		else if( token == "-o" )
		{
			if(c+1 < argc)
				outputFilename = argv[++c];
		}
		else if( token == "-c" )
		{
			useColor = true;
			if(c+1 < argc)
				color = *argv[++c];
		}
		else if( token == "-f" )
		{
			if(c+1 < argc)
				filterName = argv[++c];

			bool done2 = false;
			while( !done2 )
			{
				int nextToken2 = c+1;
				if(nextToken2 >= argc)
					done2 = true;
				else
				{
					const char * token2 = argv[nextToken2];
					if( strpbrk (token2, "-ifpoc") != NULL )
						done2 = true;
					else
					{
						params.push_back( atof(token2) );
						++c;
					}
				}
			}
		}

		if(c+1 >= argc)
			done = true;
		else
		{
			++c;
			token = argv[c];
		}
	}

	ftk::Preprocess * prep;

	if(useColor)
	{
		//Load Image
		typedef itk::ImageFileReader< ftk::Preprocess::RGBImageType3D >  ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( inputFilename );
		try
		{
			reader->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "ITK Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
			return EXIT_FAILURE;
		}

		prep = new ftk::Preprocess( reader->GetOutput(), color );
	}
	else
	{
		//Load Image
		typedef itk::ImageFileReader< ftk::Preprocess::ImageType3D >  ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( inputFilename );
		try
		{
			reader->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "ITK Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
			return EXIT_FAILURE;
		}

		prep = new ftk::Preprocess( reader->GetOutput() );
	}


	if(pipeName != "")
	{
		pipeline(pipeName, prep);
	}
	else if( filterName != "" )
	{
	}

	//Save Image:
	typedef itk::ImageFileWriter< ftk::Preprocess::ImageType3D > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( outputFilename );
	writer->SetInput( prep->GetImage() );
	try
    {
		writer->Update();
    }
	catch( itk::ExceptionObject & err ) 
    { 
		std::cout << "ITK Exception caught !" << std::endl; 
		std::cout << err << std::endl; 
		return EXIT_FAILURE;
    }

	delete prep;

	std::cerr << "PRESS ENTER TO EXIT\n";
	getchar();

	return EXIT_SUCCESS;
}

void pipeline(std::string pipeName, ftk::Preprocess * prep)
{
	TiXmlDocument doc;
	if ( !doc.LoadFile( pipeName.c_str() ) )
		return;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Preprocess" ) != 0 )
		return;

	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "LaplacianOfGaussian" ) == 0 )
		{
			int sigma = 10;
			parentElement->QueryIntAttribute("sigma",&sigma);
			std::cout << "Starting LOG...";
			prep->LaplacianOfGaussian(sigma);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "InvertIntensity" ) == 0 )
		{
			std::cout << "Starting InvertIntensity...";
			prep->InvertIntensity();
			std::cout << "done\n";
		}
		else if( strcmp( parent, "OtsuBinarize" ) == 0 )
		{
			int numThresh = 3, numFore = 2;
			parentElement->QueryIntAttribute("num_thresholds",&numThresh);
			parentElement->QueryIntAttribute("num_in_foreground",&numFore);
			std::cout << "Starting OtsuBinarize...";
			prep->OtsuBinarize(numThresh,numFore);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "RemoveConnectedComponents" ) == 0 )
		{
			int minObjSize = 1000;
			parentElement->QueryIntAttribute("minObjSize", &minObjSize);
			std::cout << "Starting RemoveCC...";
			prep->RemoveConnectedComponents(minObjSize);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "BinaryThinning" ) == 0 )
		{
			std::cout << "Starting BinaryThinning...";
			prep->BinaryThinning();
			std::cout << "done\n";
		}
		else if( strcmp( parent, "DanielssonDistanceMap" ) == 0 )
		{
			std::cout << "Starting DanielssonDistanceMap...";
			prep->DanielssonDistanceMap();
			std::cout << "done\n";
		}
		else if( strcmp( parent, "MedianFilter" ) == 0 )
		{
			int radiusX=2, radiusY=3, radiusZ=0;
			parentElement->QueryIntAttribute("radiusX",&radiusX);
			parentElement->QueryIntAttribute("radiusY",&radiusY);
			parentElement->QueryIntAttribute("radiusZ",&radiusZ);
			std::cout << "Starting MedianFilter...";
			prep->MedianFilter(radiusX,radiusY,radiusZ);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "MinErrorThresholding" ) == 0 )
		{
			float alpha_B, alpha_F, P_I;
			std::cout << "Starting MinErrorThresholding...";
			prep->MinErrorThresholding(&alpha_B, &alpha_F, &P_I);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "GraphCutBinarize" ) == 0 )
		{
			std::cout << "Starting GraphCutBinarize...";
			prep->GraphCutBinarize();
			std::cout << "done\n";
		}
		else if( strcmp( parent, "OpeningFilter" ) == 0 )
		{
			int radius=3;
			parentElement->QueryIntAttribute("radius",&radius);
			std::cout << "Starting OpeningFilter...";
			prep->OpeningFilter(radius);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "ClosingFilter" ) == 0 )
		{
			int radius=3;
			parentElement->QueryIntAttribute("radius",&radius);
			std::cout << "Starting ClosingFilter...";
			prep->ClosingFilter(radius);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "CannyEdgeDetection" ) == 0 )
		{
			float variance=2.0;
			float upperThreshold = 6;
			float lowerThreshold = 3;
			parentElement->QueryFloatAttribute("variance",&variance);
			parentElement->QueryFloatAttribute("upperThreshold",&upperThreshold);
			parentElement->QueryFloatAttribute("lowerThreshold",&lowerThreshold);
			std::cout << "Starting CannyEdgeDetection...";
			prep->CannyEdgeDetection(variance, upperThreshold, lowerThreshold);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "DiscreteGaussian") == 0 )
		{
			float varX=1.0, varY=1.0, varZ=1.0, maxError=0.1;
			parentElement->QueryFloatAttribute("varX",&varX);
			parentElement->QueryFloatAttribute("varY",&varY);
			parentElement->QueryFloatAttribute("varZ",&varZ);
			parentElement->QueryFloatAttribute("maxError",&maxError);
			std::cout << "Starting DiscreteGaussianFilter...";
			prep->DiscreteGaussianFilter(varX, varY, varZ, maxError);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "SaveVTKPoints") == 0 )
		{
			const char * filename = parentElement->Attribute("filename");
			int min=255, max=255;
			parentElement->QueryIntAttribute("min",&min);
			parentElement->QueryIntAttribute("max",&max);
			std::cout << "Saving VTK Points...";
			if(!filename)
				prep->SaveVTKPoints("points.vtk", min, max);
			else
				prep->SaveVTKPoints(filename, min, max);
			std::cout << "done\n";
		}

		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();
}

void usage(const char *funcName)
{
	std::cout << "USAGE:\n";
	std::cout << " " << funcName << " -i InputImage\n";
	//std::cout << "  -f FunctionName param1 param2 param3 ... paramN"
	std::cout << "  -o outputFile -p pipeFile -c color\n";
	std::cout << "    Default outputFile is prep_out.tif\n";
	std::cout << "    Default assumption is grayscale input, unless -c is used:\n";
	std::cout << "    -c can be 'r','g','b',or'w'\n";
}