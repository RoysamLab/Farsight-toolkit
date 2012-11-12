/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
#include "ftkProjectProcessor.h"

void usage(const char *funcName);

int main(int argc, char *argv[])
{ 
	if( argc < 5 )
	{
		usage(argv[0]);
		std::cerr << "PRESS ENTER TO EXIT\n";
		getchar();
		return EXIT_FAILURE;
	}

	std::string MyName        = argv[0];					//In case the CWD is not the path of the executable
	std::string inputFilename = argv[1];					// Name of the input image;
	std::string labelFilename = argv[2];					// Name of the label image to apply
	std::string tableFilename = argv[3];					// Name of the table file;
	std::string definitionFilename = argv[4];				// Name of the process definition file
	int numThreads = 0;
	if( argc == 6 )
		numThreads = atoi(argv[5]);					// Number of threads to run tasks built with openmp

	//Try to load the input image:
	ftk::Image::Pointer myImg = NULL;
	if( ftk::GetExtension(inputFilename) == "xml" )
	{			
		myImg = ftk::LoadXMLImage(inputFilename);
	}
	else
	{
		myImg = ftk::Image::New();
		if( !myImg->LoadFile(inputFilename) )
		{
			std::cout<<"Could not load input image\n";
			myImg = NULL;
		}
	}
	
	if(!myImg)
	{
		std::cerr << "COULD NOT LOAD INPUT IMAGE!!\n";
		usage(argv[0]);
		std::cerr << "PRESS ENTER TO EXIT\n";
		getchar();
		return EXIT_FAILURE;
	}
	
	//Try to load the Label image:
	ftk::Image::Pointer labImg;// = NULL;
	if( ftk::FileExists(labelFilename) )
	{
		if( ftk::GetExtension(labelFilename) == "xml" )
		{
			labImg = ftk::LoadXMLImage(labelFilename);
		}
		else
		{
			labImg = ftk::Image::New();
			if( !labImg->LoadFile(labelFilename) )
			{
				std::cout<<"Could not load label image\n";
				labImg = NULL;
			}
		}
	}

	//Try to load the table:
	vtkSmartPointer<vtkTable> table = NULL;
	if( ftk::FileExists(tableFilename) )
	{
		table = ftk::LoadTable(tableFilename);
	}

	//Load up the definition
	ftk::ProjectDefinition projectDef;
	if( !projectDef.Load(definitionFilename) )
	{
		std::cerr << "COULD NOT LOAD PROCESS DEFINITION FILE!!\n";
		usage(argv[0]);
		std::cerr << "PRESS ENTER TO EXIT\n";
		getchar();
		return EXIT_FAILURE;
	}

	//Do processing:
	ftk::ProjectProcessor * pProc = new ftk::ProjectProcessor();
	if( numThreads ) pProc->SetNumThreads( numThreads );
	pProc->SetExecPath( ftk::GetFilePath( MyName ) );
	std::cout<<"The executable says my path is: "<<ftk::GetFilePath( MyName )<<std::endl;
	pProc->SetInputImage(myImg);
	pProc->SetPath( ftk::GetFilePath(inputFilename) );
	if(labImg)
		pProc->SetOutputImage(labImg);
	if(table)
		pProc->SetTable(table);
	pProc->SetDefinition(&projectDef);
	pProc->Initialize();

	while(!pProc->DoneProcessing())
		pProc->ProcessNext();

	labImg = pProc->GetOutputImage();
	table = pProc->GetTable();

	typedef ftk::ProjectProcessor::LabelImageType::Pointer LabelImagePointer;
	std::map< std::string, LabelImagePointer > myClassImageMap = pProc->GetClassImageMap();
	std::map< std::string, vtkSmartPointer<vtkTable> > myClassCentroidMap = pProc->GetClassCentroidMap();
	std::string myFilename = inputFilename;
	unsigned extension = ftk::GetExtension(inputFilename).size()+1;
	std::string::iterator it;
	it = myFilename.end() - extension;
	myFilename.erase(it, it+extension);

	//Save results:
	if( ftk::GetExtension(labelFilename) == "xml" )
		ftk::SaveXMLImage(labelFilename, labImg);
	else
		labImg->SaveChannelAs(0, ftk::SetExtension(labelFilename, ""), ftk::GetExtension(labelFilename));

	ftk::SaveTable( tableFilename, table );

	typedef itk::ImageFileWriter< ftk::ProjectProcessor::LabelImageType > LabelWriterType;
	std::map< std::string, ftk::ProjectProcessor::LabelImageType::Pointer >::iterator classImageMapIter;
	for(classImageMapIter = myClassImageMap.begin(); classImageMapIter != myClassImageMap.end(); ++classImageMapIter)
	{
		std::string className = classImageMapIter->first;
		std::string classImageFileName = myFilename + "_" + className + ".nrrd";
		LabelImagePointer classImage = classImageMapIter->second;
		LabelWriterType::Pointer writer = LabelWriterType::New();
		writer->SetFileName(classImageFileName);
		writer->SetInput(classImage);
		writer->Update();
	}

	std::map< std::string, vtkSmartPointer<vtkTable> >::iterator classCentroidMapIter;
	for(classCentroidMapIter = myClassCentroidMap.begin(); classCentroidMapIter != myClassCentroidMap.end(); ++classCentroidMapIter)
	{
		std::string className = classCentroidMapIter->first;
		std::string classCentroidsFileName = myFilename + "_" + className + "_centroids.txt";
		vtkSmartPointer<vtkTable> centroid_table = classCentroidMapIter->second;

		ofstream outFile; 
		outFile.open(classCentroidsFileName.c_str(), ios::out | ios::trunc );
		if ( !outFile.is_open() )
		{
			std::cerr << "Failed to Load Document: " << outFile << std::endl;
			return false;
		}
		//Write out the features:
		for(int row = 0; row < (int)centroid_table->GetNumberOfRows(); ++row)
		{
			outFile << centroid_table->GetValue(row,0).ToInt() << "\t" ;
			outFile << centroid_table->GetValue(row,1).ToInt() << "\t" ;
			outFile << centroid_table->GetValue(row,2).ToInt() << "\t" ;
			outFile << "\n";			
		}
		outFile.close();
	}

	projectDef.Write(definitionFilename);
	
	delete pProc;

	return EXIT_SUCCESS;

}

void usage(const char *funcName)
{
	std::cout << "USAGE:\n";
	std::cout << " " << funcName << " InputImage LabelImage Table ProcessDefinition (Optional)NumThreads\n";
	std::cout << "  First four inputs are filenames\n";
}
