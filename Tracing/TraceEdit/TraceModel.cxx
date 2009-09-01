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

#include <vector>
#include <string>

#include <QtGui/QStandardItemModel>
#include <QtGui/QItemSelectionModel>

#include "TraceLine.h"
#include "TraceModel.h"

TraceModel::TraceModel(std::vector<TraceLine*> trace_lines, std::vector<std::string> FeatureHeaders)
{  
	this->Model = new QStandardItemModel(0, 0, this);
	this->SelectionModel = new QItemSelectionModel(this->Model, this);
//standard headers	
	this->headers.push_back("ID");
	this->headers.push_back("Trace Size");
	this->headers.push_back("Type");
	this->headers.push_back("Parent");
	if (FeatureHeaders.size() >=1)
	{
		for (unsigned int i = 0; i< FeatureHeaders.size(); i++)
		{
			this->headers.push_back(FeatureHeaders[i].c_str());
		}
	}
	this->SetTraces(trace_lines);
}
void TraceModel::SetTraces(std::vector<TraceLine*> trace_lines)
{
	this->TraceLines = trace_lines;
	this->SyncModel();
}
void TraceModel::SetupHeaders()
{	
	int numHeaders = this->headers.size();
	this->Model->setColumnCount(numHeaders);
	for(int i=0; i<(int)headers.size(); ++i)
    {
		this->Model->setHeaderData(i, Qt::Horizontal, 
			this->headers.at(i));
    }
}

void TraceModel::SyncModel()
{
	if(this->GetTraces().size() == 0)
	{
		return;
	}
}
