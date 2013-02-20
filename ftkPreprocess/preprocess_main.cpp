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

void usage(const char *funcName);

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
	std::string outputFilename = "";				// Name of the output image;
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
		prep->RunPipe(pipeName);
	}
	else if( filterName != "" )
	{
	}

	if( outputFilename != "" )
	{
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
		}
	}

	delete prep;

	//std::cerr << "PRESS ENTER TO EXIT\n";
	//getchar();

	return EXIT_SUCCESS;
}

void usage(const char *funcName)
{
	std::cout << "USAGE:\n";
	std::cout << " " << funcName << " -i InputImage\n";
	//std::cout << "  -f FunctionName param1 param2 param3 ... paramN"
	std::cout << "  -o outputFile -p pipeFile -c color\n";
	std::cout << "    If -o is not used the output image will not be saved\n";
	std::cout << "    Default assumption is grayscale input, unless -c is used:\n";
	std::cout << "    -c can be 'r','g','b',or'w'\n";
}