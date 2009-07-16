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
  int numRows = gaps.size();
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
  std::vector< std::vector< double > > data;
  std::vector<TraceGap*>::iterator gapItr = gaps.begin();
  for(gapItr = gaps.begin(); gapItr != gaps.end(); ++gapItr)
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
  this->setColumnCount(headers.size());

  for(int i=0; i<(int)headers.size(); ++i)
    {
    this->setHeaderData(i, Qt::Horizontal, headers.at(i));
    }

  for (int row=0; row<(int)numRows; ++row)
    {
    this->insertRow(row);
    for(int col=0; col<(int)headers.size(); ++col)
      {
      this->setData(this->index(row, col), data.at(row).at(col));
      }
    }
}

////////////////////////////////////////////////////////////////////////////////
MergeModel::~MergeModel()
{
}

////////////////////////////////////////////////////////////////////////////////
QStandardItemModel* MergeModel::GetModel()
{
  return this->Model;
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
