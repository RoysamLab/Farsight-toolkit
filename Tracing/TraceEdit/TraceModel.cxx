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

#include "TraceLine.h"
#include "TraceModel.h"

TraceModel::TraceModel(std::vector<TraceLine*> trace_lines, std::vector<std::string> FeatureHeaders)
{  
	this->DataTable = vtkSmartPointer<vtkTable>::New();
	this->Selection = new ObjectSelection();
	this->TraceClusterSelection = new SelectiveClustering();
	//this->TraceClusterManager = new ClusterManager();
	//this->TraceClusterManager->setClusteringModel(this->TraceClusterSelection );
	//this->TraceClusterManager->setObjectSelection(this->Selection);
	//this->TraceClusterManager->setManagerTitle("Trace Cluster Manager");

	//standard headers	
	this->stdHeaders();
	this->additionalFeatureHeaders.clear();
	//RPI xml currently unused
	/*if (FeatureHeaders.size() >=1)
	{
	for (int i = 0; i< (int)FeatureHeaders.size(); i++)
	{
	this->headers.push_back(FeatureHeaders[i].c_str());
	}
	}*/
	this->NumFeatures = (int)this->headers.size();
	this->SetupHeaders();
	this->SetTraces(trace_lines);
}

TraceModel::TraceModel(std::vector<TraceLine*> trace_lines)
{	this->DataTable = vtkSmartPointer<vtkTable>::New();	
this->Selection = new ObjectSelection();
this->TraceClusterSelection = new SelectiveClustering();
//this->TraceClusterManager = new ClusterManager();
//this->TraceClusterManager->setClusteringModel(this->TraceClusterSelection );
//this->TraceClusterManager->setObjectSelection(this->Selection);
//this->TraceClusterManager->setManagerTitle("Trace Cluster Manager");
//standard headers	
this->additionalFeatureHeaders.clear();
this->stdHeaders();
this->NumFeatures = (int)this->headers.size();
this->SetTraces(trace_lines);
}

TraceModel::~TraceModel()
{
	//this->TraceClusterManager->CloseClusterObjectTables();
	delete this->Selection;
	//delete this->TraceClusterManager;
	delete this->TraceClusterSelection;
}

void TraceModel::stdHeaders()
{
	this->headers.push_back("ID");
	this->headers.push_back("Level");
	this->headers.push_back("Type");
	this->headers.push_back("Root ID");
	this->headers.push_back("Path To Root");
	this->headers.push_back("Parent");
	this->headers.push_back("D To Parent");
	this->headers.push_back("# of Bits");
	this->headers.push_back("Num Children");
	this->headers.push_back("Path Length");
	this->headers.push_back("Euclidean Length");
	this->headers.push_back("Trace Density");
	this->headers.push_back("Tortuosity");
	this->headers.push_back("Radius");
	this->headers.push_back("Section Area");
	this->headers.push_back("Volume");
	this->headers.push_back("Surface Area");
	this->headers.push_back("BurkTaper");
	this->headers.push_back("HillmanTaper");
	this->headers.push_back("Bif amp Local");
	this->headers.push_back("Bif amp Remote");
	this->headers.push_back("Bif tilt Local");
	this->headers.push_back("Bif tilt Remote");
	this->headers.push_back("Bif torque Local");
	this->headers.push_back("Bif torque Remote");
	this->headers.push_back("Azimuth");
	this->headers.push_back("Elevation");
	this->headers.push_back("Terminal Degree");
	this->headers.push_back("Is Leaf");
}

void TraceModel::AddFeatureHeader(std::string NewFeatureHeader)
{
	/*!
	*
	*/
	this->additionalFeatureHeaders.push_back(NewFeatureHeader);
	this->headers.push_back(QString(NewFeatureHeader.c_str()));
}

void TraceModel::SetTraces(std::vector<TraceLine*> trace_lines)
{
	this->TraceLines.clear();
	this->TraceLines = trace_lines;
	this->SyncModel();
}
void TraceModel::SetupHeaders()
{	
	int numHeaders = (int)this->headers.size();
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	for(int i=0; i < numHeaders; ++i)
	{		
		column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( this->headers.at(i).toStdString().c_str() );
		this->DataTable->AddColumn(column);
	}
}
ObjectSelection * TraceModel::GetObjectSelection()
{
	return this->Selection;
}
void TraceModel::SyncModel()
{
	if(this->GetTraces().size() == 0)
	{
		return;
	}
	this->DataTable->Initialize();	
	this->Selection->clear();
	this->SetupHeaders();
	this->TraceIDLookupMAP.clear();
	for (int i = 0; i < (int)this->TraceLines.size(); ++i)
	{
		vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();

		DataRow->InsertNextValue(this->TraceLines.at(i)->GetId());
		DataRow->InsertNextValue((int)this->TraceLines.at(i)->GetType());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetLevel());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetRootID());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetPathLength());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetParentID(0));
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetDistToParent());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetSize());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetBranchPointer()->size());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetLength());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetEuclideanLength());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetBitDensity());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetFragmentationSmoothness());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetRadii());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetSectionArea());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetVolume());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetSurfaceArea());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetBurkTaper());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetHillmanTaper());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetBifAmplLocal());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetBifAmplRemote());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetBifTiltLocalAvg());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetBifTiltRemoteAvg());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetBifTorqueLocalAvg());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetBifTorqueRemoteAvg());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetAzimuth());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetElevation());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetTerminalDegree());
		DataRow->InsertNextValue((int)this->TraceLines.at(i)->isLeaf());
		/*for (int j = 0; j < (int)this->TraceLines.at(i)->Features.size(); ++j)
		{
		DataRow->InsertNextValue(this->TraceLines.at(i)->Features.at(j));
		}*/
		if (this->additionalFeatureHeaders.size() > 0)
		{
			for (int j = 0; j < (int)this->additionalFeatureHeaders.size();j++)
			{
				vtkVariant tempFeature = this->TraceLines.at(i)->GetTraceFeature(this->additionalFeatureHeaders[j]);
				DataRow->InsertNextValue(tempFeature);
			}
		}
		this->DataTable->InsertNextRow(DataRow);
		this->TraceIDLookupMAP[this->TraceLines.at(i)->GetId()] = this->TraceLines.at(i);
	}//end for traces.size  
	//this->MapTracesToRows();
	this->TraceClusterSelection->SetObjectTable(this->DataTable);
	//this->TraceClusterManager->setVisible(true);
}
vtkSmartPointer<vtkTable> TraceModel::getDataTable()
{
	return this->DataTable;
}
void TraceModel::SelectByIDs(int ID)
{
	this->Selection->toggle((long) ID);
	emit selectionChanged();
}
void TraceModel::SelectByIDs(std::set<long int> ID)
{
	this->Selection->add(ID);
	emit selectionChanged();
}
void TraceModel::SelectByIDs(std::vector<int> IDs)
{
	this->blockSignals(1);
	std::set<long int> ID;
	for (unsigned int i = 0; i < IDs.size(); i++)
	{
		ID.insert((long)IDs.at(i));
	}
	this->Selection->add(ID);
	this->blockSignals(0);
	emit selectionChanged();
}
void TraceModel::SetSelectionByIDs(std::vector<int> IDs)
{
	this->blockSignals(1);
	std::set<long int> ID;
	for (unsigned int i = 0; i < IDs.size(); i++)
	{
		ID.insert((long)IDs.at(i));
	}
	this->Selection->add(ID);
	this->blockSignals(0);
}
void TraceModel::SetSelectionByIDs(std::set<long int> ID)
{
	this->Selection->select(ID);
}
std::vector<int> TraceModel::GetSelectedIDs()
{
	std::vector<int> SelectedIDs;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{		
		SelectedIDs.push_back(*it);
	}
	return SelectedIDs;
}
std::vector<TraceLine*> TraceModel::GetSelectedTraces()
{
	std::vector<TraceLine*> selectedTrace;
	std::vector<int> IDList = this->GetSelectedIDs();

	//For printing out all the traces selected to the console
	//std::cout << "TraceModel IDs selected: ";
	//for (int i = 0; i < IDList.size(); i++)
	//	std::cout << IDList[i] << " ";
	//std::cout << std::endl;

	//Search for traces
	for ( unsigned int i = 0; i< IDList.size(); i++)
	{
		this->TraceIDLookupIter = this->TraceIDLookupMAP.find(IDList[i]);
		if (this->TraceIDLookupIter != this->TraceIDLookupMAP.end())
		{
			selectedTrace.push_back((*this->TraceIDLookupIter).second);
		}
	}//finished with id search
	return selectedTrace;
}
std::vector<TraceLine*> TraceModel::GetSelectedRoots()
{
	std::vector<TraceLine*> roots;
	std::vector<int> IDList;
	std::vector<TraceLine*> selectedTrace = this->GetSelectedTraces();
	for (unsigned int i = 0; i < selectedTrace.size(); i++)
	{
		int newRoot = selectedTrace.at(i)->GetRootID();
		bool found = false;
		unsigned int j= 0;
		while( !found && (j < IDList.size()))
		{
			if (IDList.at(j) == newRoot)
			{
				found = true;
			}
			else
			{
				j++;
			}
		}
		if (!found)
		{
			IDList.push_back(newRoot);
		}
	}
	for ( unsigned int i = 0; i< IDList.size(); i++)
	{
		this->TraceIDLookupIter = this->TraceIDLookupMAP.find(IDList[i]);
		if (this->TraceIDLookupIter != this->TraceIDLookupMAP.end())
		{
			roots.push_back((*this->TraceIDLookupIter).second);
		}
	}//finished with id search
	return roots;
}

void TraceModel::CloseClusterManager()
{
	//this->TraceClusterManager->close();
}

