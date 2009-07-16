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

//#include <QtGui>
//#include "TraceObject.h"
#include "TraceGap.h"
#include "TraceLine.h"
#include "MergeModel.h"

////////////////////////////////////////////////////////////////////////////////
MergeModel::MergeModel(std::vector<TraceGap*> gaps)
{
  this->SetupHeaders();
  this->SetTraceGaps(gaps);
  this->SyncModel();
}

////////////////////////////////////////////////////////////////////////////////
void MergeModel::SyncModel()
{
  if(this->GetTraceGaps().size() == 0)
    {
    return;
    }

  std::vector< std::vector< double > > data;
  std::vector<TraceGap*>::iterator gapItr = this->GetTraceGaps().begin();
  for(gapItr = this->GetTraceGaps().begin();
      gapItr != this->GetTraceGaps().end();
      ++gapItr)
    {
    std::vector<double> row;
    row.push_back((*gapItr)->compID);
    row.push_back((*gapItr)->Trace1->GetId());
    row.push_back((*gapItr)->Trace2->GetId());
    row.push_back((*gapItr)->dist);
    row.push_back((*gapItr)->angle);
    row.push_back((*gapItr)->maxdist);
    row.push_back((*gapItr)->length);
    row.push_back((*gapItr)->smoothness);
    row.push_back((*gapItr)->cost);
    data.push_back(row);
    }


  for (int row=0; row<(int)this->GetTraceGaps().size(); ++row)
    {
    //create a new row
    this->insertRow(row);
    //insert the data for a gap in this row
    for(int col=0; col < this->columnCount(); ++col)
      {
      this->setData(this->index(row, col), data.at(row).at(col));
      }
    //map gap id to row number
    }
}

////////////////////////////////////////////////////////////////////////////////
MergeModel::~MergeModel()
{
}

////////////////////////////////////////////////////////////////////////////////
void MergeModel::SetupHeaders()
{
  std::vector<QString> headers;
  headers.push_back("ID");
  headers.push_back("Trace 1");
  headers.push_back("Trace 2");
  headers.push_back("Gap");
  headers.push_back("Angle");
  headers.push_back("D");
  headers.push_back("Length");
  headers.push_back("Smoothness");
  headers.push_back("Cost");
  int numHeaders = headers.size();
  for(int i=0; i<(int)headers.size(); ++i)
    {
    this->setHeaderData(i, Qt::Horizontal, headers.at(i));
    }
  this->setColumnCount(numHeaders);
}

////////////////////////////////////////////////////////////////////////////////
void MergeModel::SetTraceGaps(std::vector<TraceGap *> gaps)
{
  this->TraceGaps = gaps;
}

////////////////////////////////////////////////////////////////////////////////
std::vector<TraceGap *> MergeModel::GetTraceGaps()
{
  return this->TraceGaps;
}

////////////////////////////////////////////////////////////////////////////////
int MergeModel::GetNumFeatures()
{ 
  return this->NumFeatures; 
}

////////////////////////////////////////////////////////////////////////////////
int MergeModel::GetNumGaps()
{ 
  return this->NumGaps; 
}

////////////////////////////////////////////////////////////////////////////////
void MergeModel::deleteTrigger(void)
{
}

////////////////////////////////////////////////////////////////////////////////
void MergeModel::mergeTrigger(void)
{
}
