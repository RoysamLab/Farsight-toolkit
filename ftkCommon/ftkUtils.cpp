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

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include "ftkUtils.h"

namespace ftk
{

bool FileExists(std::string filename)
{
	FILE * pFile = fopen (filename.c_str(),"r");
	if (pFile==NULL)
		return false;
	fclose (pFile);
	return true;
}

//Add new line to the file with the given text
bool AppendTextFile(std::string filename, std::string text)
{
	ofstream outFile; 
	outFile.open(filename.c_str(), ios::app);
	if ( !outFile.is_open() )
		return false;

	outFile << text << "\n";

	outFile.close();
	return true;
}

std::string NumToString(double d)
{
	std::stringstream out;
	out << std::setprecision(2) << std::fixed << d;	//Default is to use 2 decimal places
	return out.str();
}

std::string NumToString(int i)
{
	std::stringstream out;
	out << i ;	 
	return out.str();
}

std::string NumToString(double d, int p)
{
	std::stringstream out;
	out << std::setprecision(p) << std::fixed << d;	
	return out.str();
}

std::string TimeStamp()
{
	time_t rawtime;
	struct tm *timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	std::string dt = asctime(timeinfo);
	size_t end = dt.find('\n');
	dt.erase(end);
	return dt;
}

bool SaveTable(std::string filename, vtkSmartPointer<vtkTable> table)
{
	if(!table)
		return false;

	if(filename == "")
		return false;

	//This function writes the features to a text file
	ofstream outFile; 
	outFile.open(filename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	//Write the headers:
	for(int c=0; c<table->GetNumberOfColumns(); ++c)
	{
		outFile << table->GetColumnName(c) << "\t";
	}
	outFile << "\n";
	//Write out the features:
	for(int row = 0; row < table->GetNumberOfRows(); ++row)
	{
		for(int c=0; c < table->GetNumberOfColumns(); ++c)
		{
			outFile << ftk::NumToString( table->GetValue(row,c).ToFloat() ) << "\t";
		}
		outFile << "\n";
	}
	outFile.close();
	return true;
}

vtkSmartPointer<vtkTable> LoadTable(std::string filename)
{
	if( !FileExists(filename.c_str()) )
		return NULL;

	const int MAXLINESIZE = 1024;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	//Open the file:
	ifstream inFile; 
	inFile.open( filename.c_str() );
	if ( !inFile.is_open() )
		return NULL;

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();	

	//LOAD THE HEADER INFO:
	inFile.getline(line, MAXLINESIZE);
	char * pch = strtok (line," \t");
	while (pch != NULL)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( pch );
		table->AddColumn(column);
		pch = strtok (NULL, " \t");
	}

	//LOAD THE DATA:
	inFile.getline(line, MAXLINESIZE);
	while ( !inFile.eof() ) //Get all values
	{
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		char * pch = strtok (line," \t");
		while (pch != NULL)
		{
			row->InsertNextValue( vtkVariant( atof(pch) ) );
			pch = strtok (NULL, " \t");
		}
		table->InsertNextRow(row);
		inFile.getline(line, MAXLINESIZE);
	}
	inFile.close();
	
	return table;
}

ftk::Image::Pointer LoadXMLImage(std::string filename)
{
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return false;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Image" ) != 0 )
		return false;

	std::vector<std::string> files;
	std::vector<std::string> chName;
	std::vector<unsigned char> color;

	//Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "file" ) == 0 )
		{
			files.push_back( parentElement->GetText() );
			chName.push_back( parentElement->Attribute("chname") );
			color.push_back( atoi(parentElement->Attribute("r")) );
			color.push_back( atoi(parentElement->Attribute("g")) );
			color.push_back( atoi(parentElement->Attribute("b")) );
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();

	ftk::Image::Pointer img = ftk::Image::New();
	if(!img->LoadFilesAsMultipleChannels(files,chName,color))	//Load for display
	{
		img = NULL;
	}
	return img;
}

bool SaveXMLImage(std::string filename, ftk::Image::Pointer image)
{
	size_t pos = filename.find_last_of("/\\");
	std::string path = filename.substr(0,pos);
	pos = filename.find_last_of(".");
	std::string pre = filename.substr(0,pos);

	std::vector<std::string> names;

	//Save each channel:
	for(int i=0; i<image->GetImageInfo()->numChannels; ++i)
	{
		std::string tag = "_" + image->GetImageInfo()->channelNames.at(i);
		names.push_back( pre + tag + ".tif" );
		image->SaveChannelAs(i, pre + tag , "tif");
	}

	TiXmlDocument doc;   
 
	TiXmlElement * root = new TiXmlElement( "Image" );  
	doc.LinkEndChild( root );  
 
	for(unsigned int i=0; i<names.size(); ++i)
	  {
		TiXmlElement * file = new TiXmlElement("file");
		file->SetAttribute("chname", image->GetImageInfo()->channelNames.at(i));
		file->SetAttribute("r", NumToString(image->GetImageInfo()->channelColors.at(i).at(0)));
		file->SetAttribute("g", NumToString(image->GetImageInfo()->channelColors.at(i).at(1)));
		file->SetAttribute("b", NumToString(image->GetImageInfo()->channelColors.at(i).at(2)));
		file->LinkEndChild( new TiXmlText( names.at(i).c_str() ) );
		root->LinkEndChild(file);
	  }
	if( doc.SaveFile( filename.c_str() ) )
		return true;
	else
		return false;
}

std::string GetExtension(std::string filename)
{
	size_t pos = filename.find_last_of(".");
	std::string ext;
	if( pos == std::string::npos )
		ext = "";
	else
		ext = filename.substr(pos+1);

	return ext;
}

std::string SetExtension(std::string filename, std::string ext)
{
	std::string rName;
	size_t pos = filename.find_last_of(".");

	if(ext == "")
	{
		if( pos == std::string::npos )
			rName = filename;
		else
			rName = filename.substr(0,pos);
	}
	else
	{
		if(pos == std::string::npos)
			rName = filename + "." + ext;
		else
		{
			std::string base = filename.substr(0,pos);
			rName = base + "." + ext;
		}
	}
	return rName;
}

std::string GetFilePath(std::string f)
{
	std::string ext;
	size_t found;
	found = f.find_last_of("/\\");
	ext = f.substr(0,found);
	return ext;
}

std::string GetFilenameFromFullPath(std::string f)
{
	std::string ext;
	size_t found;
	found = f.find_last_of("/\\");
	ext = f.substr(found + 1);
	return ext;
}

std::vector<std::string> GetColumsWithString( std::string colName, vtkSmartPointer<vtkTable> table ){
	std::vector<std::string> retVect;
	for( int i=0; i<table->GetNumberOfColumns(); ++i ){
		std::string current_column;
		current_column = table->GetColumnName(i);
		if( current_column.find(colName.c_str()) != std::string::npos ){
			retVect.push_back( current_column );
		}
	}
	return retVect;	
}

std::string GetStringInCaps( std::string in_string ){
	if( !in_string.size() )
		return NULL;
	std::string out_string;
	out_string = in_string;
	std::transform(out_string.begin(), out_string.end(),out_string.begin(), ::toupper);
	return out_string;
}

}  // end namespace ftk
