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
#include "ftkProjectManager.h"

namespace ftk
{

ProjectManager::ProjectManager(const char *filename)
{
	this->fileInfo.clear();
	//read in existing project
	this->readProject(filename);
}
ProjectManager::ProjectManager()
{
	this->fileInfo.clear();
}
void ProjectManager::readProject(const char *filename)
{
	TiXmlDocument doc(filename);
	doc.LoadFile();
	if( doc.Error() )
	  {
	  std::cerr << "Error while parsing " << filename << ": " << doc.ErrorDesc() << std::endl;
	  }
	TiXmlHandle docHandle( &doc );
	TiXmlElement* sourceFile = docHandle.FirstChild("Source").FirstChild("File").Element();
	while (sourceFile)
	{
		FileInfoManager currFile;
		double tX=0,tY=0,tZ=0;
		currFile.fileName = sourceFile->Attribute("FileName");
		currFile.fileType = sourceFile->Attribute("Type");
		//translation coords. 
		sourceFile->QueryDoubleAttribute("tX",&tX);
		sourceFile->QueryDoubleAttribute("tY",&tY);
		sourceFile->QueryDoubleAttribute("tZ",&tZ);
		currFile.tx = tX;
		currFile.ty = tY;
		currFile.tz = tZ;
		this->fileInfo.push_back(currFile);
		sourceFile = sourceFile->NextSiblingElement();
	}
}
void ProjectManager::addFile(std::string fileName, std::string fileType, double x, double y, double z)
{
	FileInfoManager newFile;
	newFile.fileName = fileName;
	newFile.fileType = fileType;
	newFile.tx = x;
	newFile.ty = y;
	newFile.tz = z;
	this->fileInfo.push_back(newFile);
}
void ProjectManager::addOutputTraceFile(unsigned int i, std::string fileName)
{
	//this copies the settings from an existing file and 
	//adds the output trace file
	if (i < this->fileInfo.size())
	{
		FileInfoManager newFile;
		newFile.fileName = fileName;
		newFile.fileType = "Trace";
		//copy translation coord. 
		newFile.tx = this->fileInfo.at(i).tx;
		newFile.ty = this->fileInfo.at(i).ty;
		newFile.tz = this->fileInfo.at(i).tz;
		this->fileInfo.push_back(newFile);
	}
}
bool ProjectManager::writeProject(const char *filename)
{
	TiXmlDocument doc;
	TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
	doc.LinkEndChild( decl );
	TiXmlElement * source = new TiXmlElement( "Source" );
	doc.LinkEndChild(source);
	for (unsigned int i = 0; i < this->fileInfo.size(); i++)
	{
		TiXmlElement * sourceFile = new TiXmlElement( "File" );
		sourceFile->SetAttribute("FileName", this->fileInfo.at(i).fileName);
		sourceFile->SetAttribute("Type", this->fileInfo.at(i).fileType);
		sourceFile->SetDoubleAttribute("tX", this->fileInfo.at(i).tx);
		sourceFile->SetDoubleAttribute("tY", this->fileInfo.at(i).ty);
		sourceFile->SetDoubleAttribute("tZ", this->fileInfo.at(i).tz);
		source->LinkEndChild(sourceFile);
	}//all files should be written to nodes
	if(doc.SaveFile(filename))
	{
		return true;
	}
	else
	{	
		return false;
	}
}
unsigned int ProjectManager::size()
{
	return (int)this->fileInfo.size();
}
std::string ProjectManager::GetFileName(int i)
{
	return this->fileInfo.at(i).fileName;
}
std::string ProjectManager::GetFileType(int i)
{
	return this->fileInfo.at(i).fileType;
}
double ProjectManager::GetTranslationX(int i)
{
	return this->fileInfo.at(i).tx;
}
double ProjectManager::GetTranslationY(int i)
{
	return this->fileInfo.at(i).ty;
}
double ProjectManager::GetTranslationZ(int i)
{
	return this->fileInfo.at(i).tz;
}

void ProjectManager::ReplaceTranslations(std::string fileName, bool zOnly)
{
	const int MAXLINESIZE = 1024;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	//Open the file:
	std::ifstream inFile; 
	inFile.open( fileName.c_str() );
	if ( !inFile.is_open() )
		return;

	//inFile.getline(line, MAXLINESIZE);
	while ( !inFile.eof() ) //Get all values
	{
		inFile.getline(line, MAXLINESIZE);
	
		std::string bName;
		std::vector<double> nT;

		//Read the Image Name
		char * pch = strtok (line, " ");
		if( pch != NULL )
		{
			std::string iName = pch;
			size_t found = iName.find_last_of("."); 
			bName = iName.substr(0,found);
			//std::cerr << "Base: " << bName << std::endl;
		}

		//Read the new translation values
		pch = strtok (NULL, " ");
		while (pch != NULL)
		{
			nT.push_back( atof(pch) );
			pch = strtok (NULL, " ");
		}

		if(nT.size() < 3)
		{
			std::cerr << "Expected 3 coordinates not " << nT.size() << std::endl;
			return;
		}

		//Look for the image in project and replace the values:
		for (unsigned int i = 0; i < this->fileInfo.size(); i++)
		{
			std::string fName = this->fileInfo.at(i).fileName;
			size_t found = fName.find(bName);
			if( found != std::string::npos )
			{
				if(!zOnly)
				{
					this->fileInfo.at(i).tx = nT.at(0);
					this->fileInfo.at(i).ty = nT.at(1);
				}
				this->fileInfo.at(i).tz = nT.at(2);
			}
		}

		//inFile.getline(line, MAXLINESIZE);
	}
	inFile.close();
}


} // end namespace ftk
