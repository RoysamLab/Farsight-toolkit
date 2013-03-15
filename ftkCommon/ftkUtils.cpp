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
#include <sstream>

namespace ftk
{

bool FileExists(std::string filename)
{
	/*! 
	* Check if file Exists in the directory named 
	*/
	FILE * pFile = fopen (filename.c_str(),"r");
	if (pFile==NULL)
		return false;
	fclose (pFile);
	return true;
}

bool AppendTextFile(std::string filename, std::string text)
{
	/*!
	* Add new line to the file with the given text
	*/
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
	/*!
	* Convert a Double to std::String
	* Default is to use 2 decimal places
	*/
	std::stringstream out;
	out << std::setprecision(2) << std::fixed << d;	//
	return out.str();
}

std::string NumToString(int i)
{	/*!
	* Convert a int to std::String
	*/
	std::stringstream out;
	out << i ;	 
	return out.str();
}

std::string NumToString(double d, int p)
{
	/*!
	* Convert a Double to std::String
	* set the precision decimal places
	*/
	std::stringstream out;
	out << std::setprecision(p) << std::fixed << d;	
	return out.str();
}

std::string TimeStamp()
{
	/*!
	* 
	*/
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
	/*!
	* Write a VTK Table to a file
	*/
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

bool SaveTableAppend(std::string filename, vtkSmartPointer<vtkTable> table, int id)
{	
	/*!
	* Write a VTK Table to a file
	*/
	if(!table)
		return false;

	if(filename == "")
		return false;

	//This function writes the features to a text file
	ofstream outFile; 
	outFile.open(filename.c_str(), ios::out | ios::app );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	//Write the headers:
	//for(int c=0; c<table->GetNumberOfColumns(); ++c)
	//{
	//	outFile << table->GetColumnName(c) << "\t";
	//}
	//outFile << "\n";
	//Write out the features:
	for(int row = 0; row < table->GetNumberOfRows(); ++row)
	{
		outFile << id<<"\t";
		for(int c=1; c < table->GetNumberOfColumns(); ++c)
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
	/*!
	*	Read a tab deliminated text file and create a vtkTable
	*/
	if( !FileExists(filename.c_str()) )
		return NULL;

	const int MAXLINESIZE = 102400;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	//Open the file:
	ifstream inFile; 
	inFile.open( filename.c_str() );
	if ( !inFile.is_open() )
		return NULL;

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();	

	//!LOAD THE HEADER INFO:
	/**Creates the vtk table header information to reference the columns by*/
	inFile.getline(line, MAXLINESIZE);
	char * pch = strtok (line," \t");
	int ct = 0;
	//int i = 0;
	//std::stringstream ss;
	while (pch != NULL)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		std::string col_name = pch;
		std::replace(col_name.begin(),col_name.end(),'(','_');
		std::replace(col_name.begin(),col_name.end(),')','_');
		std::replace(col_name.begin(),col_name.end(),',','_');
		column->SetName( col_name.c_str() );
		//ss<<i++;
		//column->SetName( ss.str().c_str());
		//ss.str("");
		table->AddColumn(column);
		ct++;
		pch = strtok (NULL, " \t");
	}
	std::cout<< "Column: "<<ct<<std::endl;

	//!LOAD THE DATA:
	/*!
	* Reads all the data into table
	* Note: reads in as float but stores as vtkVariant
	*/
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

vtkSmartPointer<vtkTable> LoadRotatedTable(std::string filename)
{
	/*!
	*	Read a tab deliminated text file and create a vtkTable
	*/
	if( !FileExists(filename.c_str()) )
		return NULL;

	const int MAXLINESIZE = 102400;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	//Open the file:
	ifstream inFile; 
	inFile.open( filename.c_str());
	if ( !inFile.is_open() )
		return NULL;

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();	
	int ct = 0;
	while ( !inFile.eof() ) //Get all values
	{
		inFile.getline(line, MAXLINESIZE);
		char * pch = strtok (line," \t");

		if( pch != NULL)
		{
			vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
			std::string col_name = pch;
			if( col_name == "#N/A" )
			{
				std::stringstream out;
				out<< ct;
				col_name += out.str();
			}

			column->SetName( col_name.c_str() );
			pch = strtok (NULL, " \t");
			while (pch != NULL)
			{
				column->InsertNextValue(vtkVariant( atof(pch)));
				pch = strtok (NULL, " \t");
			}

			table->AddColumn(column);
			ct++;

			if( ct % 100 == 0)
			{
				std::cout<< ct<<"\t";
			}
		}
		else
		{
			break;
		}	
	}

	inFile.close();
	std::cout<< std::endl<< "LoadRotatedTable: "<<table->GetNumberOfRows()<<"\t"<< table->GetNumberOfColumns()<<endl;
	return table;
}

vtkSmartPointer<vtkTable> LoadXYZTable(std::string filename)
{
	/*!
	*	Read a tab deliminated text file of xyz coords and create a vtkTable
	*/
	if( !FileExists(filename.c_str()) )
		return NULL;

	const int MAXLINESIZE = 10024;	
	char line[MAXLINESIZE];

	//Open the file:
	ifstream inFile; 
	inFile.open( filename.c_str() );
	if ( !inFile.is_open() )
		return NULL;

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();	

	//!Create THE HEADER INFO:
	/**Creates the vtk table header information to reference the columns by*/

	vtkSmartPointer<vtkDoubleArray> Xcolumn = vtkSmartPointer<vtkDoubleArray>::New();
	Xcolumn->SetName( "centroid_x");
	table->AddColumn(Xcolumn);

	vtkSmartPointer<vtkDoubleArray> Ycolumn = vtkSmartPointer<vtkDoubleArray>::New();
	Ycolumn->SetName( "centroid_y");
	table->AddColumn(Ycolumn);

	vtkSmartPointer<vtkDoubleArray> Zcolumn = vtkSmartPointer<vtkDoubleArray>::New();
	Zcolumn->SetName( "centroid_z");
	table->AddColumn(Zcolumn);
	//!LOAD THE DATA:
	/*!
	* Reads all the data into table
	* Note: reads in as float but stores as vtkVariant
	*/
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
vtkSmartPointer<vtkTable> AppendLoadTable(std::string filename, vtkSmartPointer<vtkTable> initialTable , double tx, double ty, double tz)
{	//!Loads and apends a feature table into montage space
	/*!
	* Filename is for the loading the new file
	* initialTable is the table to append to
	* tx, ty, tz are the shift value to position into montage
	*/
	vtkSmartPointer<vtkTable> outputTable;
	vtkSmartPointer<vtkTable> newTable = ftk::LoadTable(filename);
	vtkVariant maxid = 0;
	if (initialTable->GetNumberOfRows() > 0)
	{
		for(vtkIdType r = 0; r < initialTable->GetNumberOfRows() ; r++ )
		{
			vtkVariant tempID = initialTable->GetValueByName( r, "ID");
			if (tempID > maxid)
			{
				maxid = tempID;
			}
		}//end for rows
	}
	for (vtkIdType rnew = 0; rnew < newTable->GetNumberOfRows() ; rnew++ )
	{
		vtkVariant tempID = newTable->GetValueByName(rnew , "ID");
		newTable->SetValueByName(rnew, "ID", vtkVariant(tempID.ToInt() + maxid.ToInt()));

		vtkVariant centroid_X = newTable->GetValueByName(rnew , "centroid_x");
		newTable->SetValueByName(rnew , "centroid_x",  vtkVariant(centroid_X.ToDouble() + tx));

		vtkVariant centroid_Y = newTable->GetValueByName(rnew , "centroid_y");
		newTable->SetValueByName(rnew , "centroid_y",  vtkVariant(centroid_Y.ToDouble() + ty));

		vtkVariant centroid_Z = newTable->GetValueByName(rnew , "centroid_z");
		newTable->SetValueByName(rnew , "centroid_z",  vtkVariant(centroid_Z.ToDouble() + tz));
	}
	if (initialTable->GetNumberOfRows() > 0)
	{
		outputTable = ftk::AppendTables(initialTable, newTable);
	}
	else
	{
		outputTable = newTable;
	}
	return outputTable;

}
std::vector< vtkSmartPointer<vtkTable> > LoadTableSeries(std::string filename)
{
	/*!
	* Loads multiple tables into a vector to access
	* Ie. Time series data tables
	*/
	std::vector< vtkSmartPointer<vtkTable> > tableVector;
	tableVector.clear();

	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return tableVector;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Table" ) != 0 )
		return tableVector;

	//Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();

	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "file" ) == 0 )
		{
			tableVector.push_back(LoadTable(std::string(reinterpret_cast<const char*>(parentElement->GetText()))));
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();
	
	return tableVector;
}



bool SaveTableSeries(std::string filename,std::vector< vtkSmartPointer<vtkTable> >  table4DImage,std::string path)
{

	std::string seriesfilename = filename ;
	TiXmlDocument doc;    
	TiXmlElement * root = new TiXmlElement( "Image" );  
	doc.LinkEndChild( root ); 
	


	for(unsigned int i = 0;	i<table4DImage.size(); ++i)				// iterate through time
	{
		// write the xml file containing the xml file names (time file names) of other xml files representing channels
		TiXmlElement * file = new TiXmlElement("file");
		std::stringstream ss;
		ss<<i;
		std::string name = path+"table_time"+ss.str()+".txt";
		SaveTable(name,table4DImage.at(i));
		file->LinkEndChild( new TiXmlText( name ) );
		std::cout<< name<<std::endl;
		root->LinkEndChild(file);
	}
	if( doc.SaveFile( seriesfilename.c_str() ) )
		return true;
	else
		return false;
}




bool SaveTableSeriesActive(std::string filename,std::vector< vtkSmartPointer<vtkTable> >  table4DImage)
{

	std::vector< vtkSmartPointer<vtkTable> > tableVector;
	tableVector.clear();

	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return false;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Table" ) != 0 )
		return false;

	//Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();

	int counter = 0;

	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "file" ) == 0 )
		{
			SaveTable(std::string(reinterpret_cast<const char*>(parentElement->GetText())),table4DImage.at(counter) );
		}
		parentElement = parentElement->NextSiblingElement();
		counter++;
	} // end while(parentElement)
	//doc.close();
	return true;
}




vtkSmartPointer<vtkTable> AppendTables(vtkSmartPointer<vtkTable> table_initial,vtkSmartPointer<vtkTable> table_new )
{
	/*!
	* Adds on table to the end of another
	*/
  //!fill the table with values
  unsigned int counter = 0;
  for(vtkIdType r = 0; r < table_new->GetNumberOfRows() ; r++ )
    {
	  table_initial->InsertNextRow(table_new->GetRow(r));
    }

	return table_initial;
}




ftk::Image::Pointer LoadXMLImage(std::string filename)
{
	/*!
	* Load a ftk::image from xml file
	*/
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return NULL;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Image" ) != 0 )
		return NULL;

	std::vector<std::string> files;
	std::vector<std::string> chName;
	std::vector<unsigned char> color;

	//!Parents we know of: datafilename,resultfilename,object,parameter
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

ftk::Image::Pointer LoadImageSeries(std::string filename)
{
	/*!
	* Loads a time series
	*/
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return NULL;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Image" ) != 0 )
		return NULL;


	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(1); //RELEASE_CONTROL mode	

	//!Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	ftk::Image::Pointer img = ftk::Image::New();
	img = LoadXMLImage(std::string(reinterpret_cast<const char*>(parentElement->GetText())));
	parentElement = parentElement->NextSiblingElement();

	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "file" ) == 0 )
		{
			img->AppendImage(LoadXMLImage(std::string(reinterpret_cast<const char*>(parentElement->GetText()))),mode,true);
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();
	return img;
}

std::vector<std::string> GetSeriesPaths(std::string filename){
	std::vector<std::string> paths;
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return paths;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Image" ) != 0 )
		return paths;

	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement){
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "file" ) == 0 ){
			std::string file_name = std::string(reinterpret_cast<const char*>(parentElement->GetText()));
			paths.push_back( file_name );
			file_name.clear();
		}
		parentElement = parentElement->NextSiblingElement();
	}

	return paths;
}


ftk::Image::Pointer LoadImageSeriesLabels(std::string filename)
{
	/*!
	* Loads a time series label images
	*/
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return NULL;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Image" ) != 0 )
		return NULL;


	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(1); //RELEASE_CONTROL mode	

	//!Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	ftk::Image::Pointer img = ftk::Image::New();

	img->LoadFile(std::string(reinterpret_cast<const char*>(parentElement->GetText())));
	if(!img)
	{	
		img = NULL;
	}
		
	parentElement = parentElement->NextSiblingElement();

	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "file" ) == 0 )
		{
			ftk::Image::Pointer tempImg = ftk::Image::New();
			tempImg->LoadFile(std::string(reinterpret_cast<const char*>(parentElement->GetText())));
			img->AppendImage(tempImg,mode,true);
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();

	std::cout<<img->GetImageInfo()->numTSlices<<std::endl;

	return img;
}
bool SaveLabelSeries(std::string seriesfilename,ftk::Image::Pointer image, std::string path)
{
	if(seriesfilename == "")
		return false;

	if(GetExtension(seriesfilename)!="xml")
		return false;
	std::vector<std::vector<std::string> > filenames;
	//if(image->GetImageInfo()->numTSlices == 1)
	//	filenames.push_back(image->GetFilenames());
	//else
		filenames = image->GetTimeChannelFilenames();

	//seriesfilename = GetFilePath(filenames.at(0))+"\\"+seriesfilename ;
	seriesfilename = path+seriesfilename ;
	TiXmlDocument doc;    
	TiXmlElement * root = new TiXmlElement( "Image" );  
	doc.LinkEndChild( root ); 
	


	for(unsigned int i = 0;	i<filenames.size(); ++i)				// iterate through time
	{
		// write the xml file containing the xml file names (time file names) of other xml files representing channels
		TiXmlElement * file = new TiXmlElement("file");
	//	std::string path = GetFilePath(filenames.at(i));
		//std::string fullname = GetFilenameFromFullPath(filenames.at(i).at(0));		// first channel is assumed to be segmented nucleus
		//std::string fname = path+"\\"+"labeled_"+fullname;
		file->LinkEndChild( new TiXmlText( filenames.at(i).at(0).c_str() ) );
		std::cout<< filenames.at(i).at(0)<<std::endl;
		root->LinkEndChild(file);
	}
	if( doc.SaveFile( seriesfilename.c_str() ) )
		return true;
	else
		return false;
}	

bool SaveImageSeries(std::string seriesfilename, ftk::Image::Pointer image, std::string path)
{
	if(!image)
		return false;

	if(seriesfilename == "")
		return false;

	if(GetExtension(seriesfilename)!="xml")
		return false;
	
	std::vector<std::vector<std::string> > filenames;
	if(image->GetImageInfo()->numTSlices == 1)
		filenames.push_back(image->GetFilenames());
	else
		filenames = image->GetTimeChannelFilenames();	

	//seriesfilename = GetFilePath(filenames.at(0).at(0))+"\\"+seriesfilename ;
	seriesfilename = path+seriesfilename ;
	TiXmlDocument doc;    
	TiXmlElement * root = new TiXmlElement( "Image" );  
	doc.LinkEndChild( root ); 

	for(unsigned int i = 0;	i<filenames.size(); ++i)				// iterate through time
	{
		// write the xml file containing the xml file names (time file names) of other xml files representing channels
		TiXmlElement * file = new TiXmlElement("file");
		std::stringstream ss;//create a stringstream
		ss << i+1;
		std::string xmlfilename = path+"series_time_"+ss.str()+".xml";
		file->LinkEndChild( new TiXmlText( xmlfilename.c_str() ) );
		root->LinkEndChild(file);

		// write the xml file defining the channels for this time point 
		TiXmlDocument subdoc;
		TiXmlElement * subroot = new TiXmlElement( "Image" );  
		subdoc.LinkEndChild( subroot ); 

		for(unsigned int j=0; j<filenames.at(i).size(); ++j)
		  {
			TiXmlElement * subfile = new TiXmlElement("file");
			subfile->SetAttribute("chname", image->GetImageInfo()->channelNames.at(j));
			subfile->SetAttribute("r", NumToString(image->GetImageInfo()->channelColors.at(j).at(0)));
			subfile->SetAttribute("g", NumToString(image->GetImageInfo()->channelColors.at(j).at(1)));
			subfile->SetAttribute("b", NumToString(image->GetImageInfo()->channelColors.at(j).at(2)));
			subfile->LinkEndChild( new TiXmlText( filenames.at(i).at(j).c_str() ) );
			subroot->LinkEndChild(subfile);
		  }
		subdoc.SaveFile( xmlfilename.c_str() );

	}
	if( doc.SaveFile( seriesfilename.c_str() ) )
		return true;
	else
		return false;
}	





bool SaveXMLImage(std::string filename, ftk::Image::Pointer image)
{
	/*!
	* save an ftk::image to a multichannel xml file
	*/
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



std::vector<Channel> ReadChannels(TiXmlElement * inputElement)
{
	std::vector<Channel> returnVector;

	TiXmlElement * channelElement = inputElement->FirstChildElement();
	while (channelElement)
	{
		const char * parent = channelElement->Value();
		if ( strcmp( parent, "channel" ) == 0 )
		{
			Channel channel;
			channel.number = atoi(channelElement->Attribute("number"));
			channel.name = channelElement->Attribute("name");
			channel.type = channelElement->Attribute("type");
			returnVector.push_back(channel);
		}
		channelElement = channelElement->NextSiblingElement();
	} // end while(channelElement)
	return returnVector;
}



//std::vector<ftk::AssociationRule> ReadAssociationRules(TiXmlElement * inputElement)
//{
//	std::vector<ftk::AssociationRule> returnVector;
//
//	TiXmlElement * parameterElement = inputElement->FirstChildElement();
//	while (parameterElement)
//	{
//		const char * parameter = parameterElement ->Value();
//		ftk::AssociationRule assocRule("");
//		if ( strcmp(parameter,"AssociationRule") == 0 )
//		{
//			assocRule.SetRuleName(parameterElement->Attribute("Name"));
//			assocRule.SetSegmentationFileNmae(parameterElement->Attribute("SegmentationSource"));
//			assocRule.SetTargetFileNmae(parameterElement->Attribute("Target_Image"));
//			assocRule.SetOutDistance(atoi(parameterElement->Attribute("Outside_Distance")));
//			assocRule.SetInDistance(atoi(parameterElement->Attribute("Inside_Distance")));
//
//			if(strcmp(parameterElement->Attribute("Use_Whole_Object"),"True")==0)
//				assocRule.SetUseWholeObject(true);
//			else
//				assocRule.SetUseWholeObject(false);
//
//			if(strcmp(parameterElement->Attribute("Use_Background_Subtraction"),"True")==0)
//				assocRule.SetUseBackgroundSubtraction(true);
//			else
//				assocRule.SetUseBackgroundSubtraction(false);
//
//			if(strcmp(parameterElement->Attribute("Use_MultiLevel_Thresholding"),"True")==0){
//				assocRule.SetUseMultiLevelThresholding(true);
//				assocRule.SetNumberOfThresholds(atoi(parameterElement->Attribute("Number_Of_Thresholds")));
//				assocRule.SetNumberIncludedInForeground(atoi(parameterElement->Attribute("Number_Included_In_Foreground")));
//			}
//			else{
//				assocRule.SetUseMultiLevelThresholding(false);
//				assocRule.SetNumberOfThresholds(1);
//				assocRule.SetNumberIncludedInForeground(1);
//			}
//
//			if(strcmp(parameterElement->Attribute("Association_Type"),"MIN")==0)
//				assocRule.SetAssocType(ASSOC_MIN);
//			else if(strcmp(parameterElement->Attribute("Association_Type"),"MAX")==0)
//				assocRule.SetAssocType(ASSOC_MAX);
//			else if(strcmp(parameterElement->Attribute("Association_Type"),"TOTAL")==0)
//				assocRule.SetAssocType(ASSOC_TOTAL);
//			else if(strcmp(parameterElement->Attribute("Association_Type"),"SURROUNDEDNESS")==0)
//				assocRule.SetAssocType(ASSOC_SURROUNDEDNESS);
//			else
//				assocRule.SetAssocType(ASSOC_AVERAGE);
//		}
//		returnVector.push_back(assocRule);
//		parameterElement = parameterElement->NextSiblingElement();
//	}
//	return returnVector;
//}




std::string GetExtension(std::string filename)
{
	/*!
	* Reads File extensions 
	*/
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
	//! Set file to specified extension
	/*!
	* Remove old extension
	* Append a specified file extension to filename
	*/
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
	/*!
	* returns the file path portion of a string
	*/
	std::string ext;
	size_t found;
	found = f.find_last_of("/\\");
	if( found != std::string::npos )
		ext = f.substr(0,found);
	else
		ext = "./";
	return ext;
}

std::string GetFilenameFromFullPath(std::string f)
{
	/*!
	* returns the file name portion of a string
	*/
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


vtkSmartPointer<vtkTable> CopyTable(vtkSmartPointer<vtkTable> featureTable )
{
	/*!
	* Make a dupilcate vtkTable 
	* "deep copy"
	*/
	vtkSmartPointer<vtkTable> table_validation  = vtkSmartPointer<vtkTable>::New();
	table_validation->Initialize();

	for(int col=0; col<featureTable->GetNumberOfColumns(); ++col)
	{	
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName(featureTable->GetColumnName(col));
		table_validation->AddColumn(column);	
	}

	for(int row = 0; row < (int)featureTable->GetNumberOfRows(); ++row)
	{	
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int c =0;c<(int)table_validation->GetNumberOfColumns();++c)
			model_data1->InsertNextValue(featureTable->GetValueByName(row,table_validation->GetColumnName(c)));
		table_validation->InsertNextRow(model_data1);
	}

	return table_validation;
}

double GetMean(std::vector<double> data)
{
	if(data.empty())
		return 0.0;

	double sum = 0;
	for(int i=0; i<(int)data.size(); ++i)
	{
		sum+=data[i];
	}
	return sum/(double)data.size();
}
double GetStd(std::vector<double> data)
{
	if(data.empty())
		return 0.0;

	double sum = 0;
	double sumx2 = 0;
	for(int i=0; i<(int)data.size(); ++i)
	{
		sum+=data[i];
		sumx2+= data[i]*data[i]; 
	}
	sum /=(double)data.size();
	sumx2 /= (double)data.size();

	return sumx2-(sum*sum);
}

}  // end namespace ftk
