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
#include "TraceProjectManager.h"

ProjectManager::ProjectManager(char *filename)
{
	this->fileInfo.clear();
	//read in existing project
	TiXmlDocument doc(filename);
	doc.LoadFile();
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
bool ProjectManager::writeProject(char *filename)
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
