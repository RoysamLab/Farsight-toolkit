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
#include <iostream>
#include <math.h>
#include <queue>
#include <set>
#include <QProgressDialog>
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkDataArray.h"
#include "vtkAppendPolyData.h"
#include "vtkSphereSource.h"
#include "vtkGlyph3D.h"
#include "vtkArrowSource.h"
#include "TraceBit.h"
#include "TraceLine.h"
#include "TraceObject.h"
#include "TraceGap.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New()

#include "tinyxml/tinyxml.h"

#define MY_ENCODING "ISO-8859-1"

TraceObject::TraceObject()
{
}

TraceObject::TraceObject(const TraceObject &T)
{
  for(unsigned int counter=0; counter<T.trace_lines.size(); counter++)
  {
    TraceLine* temp = new TraceLine();
    *temp = *(T.trace_lines[counter]);
    this->trace_lines.push_back(temp);
  }
  //this->SmallLines = T.SmallLines; FIXME: is  'SmallLines' storing some redundant information? we need to take care of it in the copy constructor.
  this->smallLineColor = T.smallLineColor;
  this->mergeLineColor = T.mergeLineColor;
}

TraceObject::~TraceObject()
{
  for(unsigned int counter=0; counter<trace_lines.size(); counter++)
    delete trace_lines[counter];
}

double TraceObject::getSmallLineColor()
{
  return this->smallLineColor;
}

double TraceObject::getMergeLineColor()
{
  return this->mergeLineColor;
}

void TraceObject::setSmallLineColor(double set)
{
  this->smallLineColor=set;
}

void TraceObject::setMergeLineColor(double set)
{
  this->mergeLineColor=set;
}

void TraceObject::Print(std::ostream &c)
{
  c<<"TraceObject:"<<std::endl;
  c<<"Size:"<<trace_lines.size()<<std::endl;
  for(unsigned int counter=0; counter<trace_lines.size(); counter++)
  {
    trace_lines[counter]->Print(c,4);
  }
}

std::vector<TraceLine*>* TraceObject::GetTraceLinesPointer()
{
  return &trace_lines;
}

std::vector<TraceLine*> TraceObject::GetTraceLines()
{
  return trace_lines;
}

bool TraceObject::ReadFromFeatureTracksFileForKymograph(char *filename,int type_offset=0)
{
  FILE * fp = fopen(filename,"r");
  if(fp == NULL)
  {
    printf("Could not open %s for reading\n",filename);
    return false;
  }
  int curr_track = -1;
  TraceLine* curr_line = NULL;
  int line_count=0;
  int curr_id = 0;
  while(!feof(fp))
  {
    int track;
    double x,y,z;
    int t;
    TraceBit tbit;
    fscanf(fp,"%d %d %lf %lf %lf",&track,&t,&x,&y,&z);//We do not represent the z axis, instead use it for coloring the tracks

    tbit.x = x; tbit.y = 2.97/2.79*y; tbit.z = t; tbit.id = z;line_count++;   
    tbit.r = 1;

    if(track!=curr_track)
    {
      curr_track = track;
      TraceLine* tline = new TraceLine();
      trace_lines.push_back(tline);
      tline->SetId(curr_id++);
      tline->SetType(1+type_offset);
      curr_line = tline;
    }
    curr_line->AddTraceBit(tbit);
  }
  printf("I read %d lines from the file %s\n",line_count,filename);
  fclose(fp);
  return true;
}

bool TraceObject::ReadFromFeatureTracksFile(char *filename,int type_offset=0)
{
  FILE * fp = fopen(filename,"r");
  if(fp==NULL)
    return false;
  int curr_track = -1;
  TraceLine* curr_line = NULL;
  int curr_id = 0;
  while(!feof(fp))
  {
    int track;
    TraceBit tbit;
    fscanf(fp,"%d %d %lf %lf %lf",&track,&tbit.id,&tbit.x,&tbit.y,&tbit.z);
    tbit.r = 1;
    if(track!=curr_track)
    {
      curr_track = track;
      TraceLine* tline = new TraceLine();
      tline->SetId(curr_id++);
      tline->SetType(1+type_offset);
      trace_lines.push_back(tline);
      curr_line = tline;
    }
    curr_line->AddTraceBit(tbit);
  }
  fclose(fp);
  return true;
}

bool TraceObject::ReadFromSWCFile(char * filename)
{

  FILE * fp = fopen(filename, "r");
  if(fp==NULL)
  {
    printf("Couldn't open file %s for parsing\n",filename);
    return false;
  }
  char buff[1024];
  
  //make an initial pass through the file to figure out how many points
  //we're dealing with
  int numPoints = 1;
  while(!feof(fp))
    {
    if(fgets(buff,1024,fp)==NULL)
      {
      break;
      }
    int pc = 0;
    while(buff[pc]==' '&&(pc<1023))
      {
      pc++;
      }
    if(buff[pc]=='#') 
      {
      continue;
      }
	int temp_id;
	sscanf(buff,"%d %*d %*f %*f %*f %*f %*d",&temp_id);
	numPoints = (numPoints>temp_id)?numPoints:temp_id;
    }
  numPoints++; // set it to 1 + maximum id in the file
  rewind(fp);

  unsigned char *child_count = (unsigned char *)malloc(numPoints * sizeof(char));
  std::set<int> criticals; // store all points who have parent = -1 or child_count[parent] > 1

  vtksys::hash_map<unsigned int,int> hash_type; // smaller hash functions only for the critical points.. saves us memory and time
  vtksys::hash_map<unsigned int,int> hash_parent;
  vtksys::hash_map<unsigned int,unsigned long long int> hash_load;

  //memset(child_count,0,sizeof(unsigned char)*100000);
  
  for(int counter=0; counter < numPoints; counter++)
    {
    child_count[counter]=0;
    }
  int id, type,parent;
  double x,y,z,r;
  int max_id = -1;
  while(!feof(fp))
    {
    if(fgets(buff,1024,fp)==NULL)
      break;
    int pc = 0;
    while(buff[pc]==' '&&(pc<1023))
      pc++;
    if(buff[pc]=='#') // ignoring comment lines for now. We have to read scale from it! TODO
      continue;

    //sscanf(buff,"%d %d %*f %*f %*f %*f %d",&id,&type,&parent);
    sscanf(buff,"%d %d %*f %*f %*f %*f %d",&id,&type,&parent);
    //printf("%d\n",id);
    if(id>max_id)//find max id
    {
      max_id=id;
    }

    if(parent >= 1)
    {
      child_count[parent]++;
    }
	hash_type[id] = type;
  }
  fclose(fp);
  //printf("I read %d lines\n",tc);
  unsigned int *child_id = (unsigned int *)malloc(numPoints * sizeof(unsigned int));
  //memset(child_id,0,sizeof(unsigned int)*(max_id));
  std::vector<TraceBit> data(max_id+1);

  fp = fopen(filename,"r");
  int tcc =0;
  while(!feof(fp))
  {
    tcc++;
    //printf("Done %d\n",tcc);
    if(fgets(buff,1024,fp)==NULL)
      break;
    int pc = 0;
    while(buff[pc]==' '&&(pc<1023))
      pc++;
    if(buff[pc]=='#') // ignoring comment lines for now. We have to read scale from it! TODO
      continue;
    sscanf(buff,"%d %d %lf %lf %lf %lf %d",&id,&type,&x,&y,&z,&r,&parent);
    TraceBit tbit;
    tbit.x=x;tbit.y=y;tbit.z=z;tbit.id=id;tbit.r =r;
    data[id] = tbit;
    if(parent!=-1)
    {
      child_id[parent] = id;
    }
    if(parent == -1)
    {
      criticals.insert(id);
      hash_parent[id] = -1;
      //printf("hash_parent[%d] = %d\n",id,hash_parent[id]);
    }
    else
    {
      if(child_count[parent]>1)
      {
        criticals.insert(id);
        hash_parent[id] = parent;
        //printf("hash_parent[%d] = %d\n", id, hash_parent[id]);
      }
	  else if(hash_type[id] != hash_type[parent])
	  {
		  criticals.insert(id);
		  hash_parent[id] = parent;
	  }
    }
  }
  fclose(fp);
  printf("about to create the data structure.. %d\n",criticals.size());
  std::set<int>::iterator iter = criticals.begin();
  int global_id_number = 1;
  while(iter != criticals.end())
  {
    TraceLine * ttemp = new TraceLine();
    ttemp->SetId(global_id_number++);
    ttemp->SetType(hash_type[*iter]);
    //ttemp->setTraceColor(1.0/ttemp->GetType());
    switch( ttemp->GetType() )
    {
    case 1:
        ttemp->setTraceColor(.75);//cyan
        break;
    case 3:
    case 4:
    case 5:
        ttemp->setTraceColor(.90);//blue
        break;
    case 7:
    case 0:
    default:
        ttemp->setTraceColor(.25); //yellow
    }
    ttemp->AddTraceBit(data[*iter]);
    int id_counter = *iter;
    while(child_count[id_counter]==1)
    {
		if(hash_type[id_counter] == hash_type[child_id[id_counter]])
		{
			id_counter = child_id[id_counter];
			ttemp->AddTraceBit(data[id_counter]);
		}
		else
			break;

    }
    hash_load[id_counter] = reinterpret_cast<unsigned long long int>(ttemp); 
    // Important: We're storing TraceLine* for points in the end of segments only.
    trace_lines.push_back(ttemp);
    iter++;
  }
  //printf("Trace_lines size = %d\n",trace_lines.size());
  
  iter = criticals.begin();
  int pc = 0;
  while(iter!= criticals.end())
  {
    //printf("trace_lines[%d] = %p\n",pc,trace_lines[pc]);
    if(hash_parent[*iter]>0)
    {
      //printf("hash_parent %d *iter %d hash_load %p\n",hash_parent[*iter],*iter,reinterpret_cast<void*>(hash_load[hash_parent[*iter]]));
      TraceLine * t = reinterpret_cast<TraceLine*>(hash_load[hash_parent[*iter]]);
      trace_lines[pc]->SetParent(t);
      
      //t->AddBranch(trace_lines[pc]);
      t->GetBranchPointer()->push_back(trace_lines[pc]);
      //if(t->GetBranchPointer()->size()>2)
      //{
      //  printf("here is the error\n");
      //}
    }
    else
    {
      trace_lines[pc]->SetParent(NULL);
    }
    pc++;
    ++iter;
  }
  std::vector<TraceLine*>::iterator cleaniter = trace_lines.end();
  cleaniter--;
  for(;cleaniter!=trace_lines.begin(); cleaniter--)
  {
    if((*cleaniter)->GetParent()!=NULL)
      cleaniter=trace_lines.erase(cleaniter);
  }
  printf("Finished loading\n");
  //Print(std::cout);
//  delete [] child_id;
  free(child_count);
  free(child_id);
  return true;
}

void TraceObject::CreatePolyDataRecursive(TraceLine* tline, vtkSmartPointer<vtkFloatArray> point_scalars, 
                      vtkSmartPointer<vtkPoints> line_points,vtkSmartPointer<vtkCellArray> line_cells)
{
  //tline->Print(std::cout);
  //scanf("%*c");
  TraceLine::TraceBitsType tbits;
  TraceLine::TraceBitsType::iterator iter=tline->GetTraceBitIteratorBegin();
  double point[3];
  unsigned int return_id;
  unsigned int cell_id;
  unsigned int old_id;

  std::vector<unsigned int>* cell_id_array=tline->GetMarkers();
  cell_id_array->clear();

  point[0] = iter->x;point[1]=iter->y;point[2]=iter->z;
  return_id = line_points->InsertNextPoint(point);
  hashp[return_id]=(unsigned long long int)tline;
  iter->marker = return_id;
  point_scalars->InsertNextTuple1(.5-tline->getTraceColor());

  //To add a line between parent line's last point and the first point in the current line
  if(tline->GetParent() != NULL)
  {
    //printf("I should not have a parent at all! why did I come here?\n");
    if(tline->GetParent()->GetTraceBitsPointer()->size()>0)
    {
      cell_id = line_cells->InsertNextCell(2);
      cell_id_array->push_back(cell_id);
      hashc[cell_id] = reinterpret_cast<unsigned long long int>(tline);
      line_cells->InsertCellPoint((--(tline->GetParent()->GetTraceBitIteratorEnd()))->marker);
      line_cells->InsertCellPoint(return_id);
    }
  }
  // Rest of the lines for the current tline 
  iter++;
  while(iter!=tline->GetTraceBitIteratorEnd())
  {
    //printf("in loop %d\n",++pc);
    old_id = return_id;
    point[0] = iter->x;point[1]=iter->y;point[2]=iter->z;
    return_id = line_points->InsertNextPoint(point);
    hashp[return_id]=(unsigned long long int)tline;
    iter->marker = return_id;

    point_scalars->InsertNextTuple1(tline->getTraceColor());
    cell_id = line_cells->InsertNextCell(2);
    cell_id_array->push_back(cell_id);
    hashc[cell_id]=reinterpret_cast<unsigned long long int>(tline);
    line_cells->InsertCellPoint(old_id);
    line_cells->InsertCellPoint(return_id);
    ++iter;
  }
  // Recursive calls to the branches if they exist 
  for(unsigned int counter=0; counter<tline->GetBranchPointer()->size(); counter++)
  {
    //printf("I should be having children too! what am I doing here?\n");
    CreatePolyDataRecursive((*tline->GetBranchPointer())[counter],point_scalars,line_points,line_cells);
  }
}

void TraceObject::CollectTraceBitsRecursive(std::vector<TraceBit> &vec,TraceLine *l)
{
  TraceLine::TraceBitsType *p = l->GetTraceBitsPointer();
  TraceLine::TraceBitsType::iterator iter = p->begin();
  TraceLine::TraceBitsType::iterator iterend = p->end();
  while(iter!=iterend)
  {
    vec.push_back(*iter);
    iter++;
  }
  std::vector<TraceLine*>* bp = l->GetBranchPointer();
  for(unsigned int counter=0; counter< bp->size(); counter++)
  {
    this->CollectTraceBitsRecursive(vec,(*bp)[counter]);
  }
}

std::vector<TraceBit> TraceObject::CollectTraceBits()
{
  std::vector<TraceBit> vec;
  for(unsigned int counter=0; counter<trace_lines.size(); counter++)
    this->CollectTraceBitsRecursive(vec,trace_lines[counter]);
  return vec;
}

struct equlli
{
  bool operator()(const unsigned long long int l1, const unsigned long long int l2) const
  {
    return l1 == l2;
  }
};
//needed because gcc doesn't have a built-in method to hash unsigned long long ints
struct hashulli
{
  size_t operator()(const unsigned long long int __x) const
  {
    return __x;
  }
  size_t operator()(const unsigned long long int __x, const unsigned long long int __y)
  {
    return __x == __y;
  }
  const static size_t bucket_size = 4;
  const static size_t min_buckets = 8;
};

bool TraceObject::WriteToSWCFile(const char *filename)
{
  FILE * fp = fopen(filename,"w");
  vtksys::hash_map<const unsigned long long int, int, hashulli> hash_dump;
  if(fp == NULL)
  {
    printf("Couldn't open %s for writing\n",filename);
    return false;
  }
  int cur_id = 1;
  std::queue<TraceLine*> q;
  for(unsigned int counter=0; counter<trace_lines.size(); counter++)
  {
    q.push(trace_lines[counter]);
  }
  while(!q.empty())
  {
    TraceLine *t = q.front();
    q.pop();
    TraceLine::TraceBitsType::iterator iter = t->GetTraceBitIteratorBegin();
    TraceLine::TraceBitsType::iterator iterend = t->GetTraceBitIteratorEnd();
    hash_dump[reinterpret_cast<unsigned long long int>(t)]=cur_id+t->GetTraceBitsPointer()->size()-1;
    if(t->GetParent()==NULL)
    {
      fprintf(fp,"%d %d %0.2lf %0.2lf %0.2lf %0.2lf %d\n",cur_id++,t->GetType(),iter->x,iter->y,iter->z,iter->r,-1);
    }
    else
    {
      fprintf(fp,"%d %d %0.2lf %0.2lf %0.2lf %0.2lf %d\n",cur_id++,t->GetType(),iter->x,iter->y,iter->z,iter->r,hash_dump[reinterpret_cast<unsigned long long int>(t->GetParent())]);
    }
    iter++;
    while(iter!=iterend)
    {
      fprintf(fp,"%d %d %0.2lf %0.2lf %0.2lf %0.2lf %d\n",cur_id,t->GetType(),iter->x,iter->y,iter->z,iter->r,cur_id-1);
      cur_id++;
      iter++;
    }
    for(unsigned int counter=0; counter<t->GetBranchPointer()->size(); counter++)
    {
      q.push((*t->GetBranchPointer())[counter]);
    }
  }
  fclose(fp);
  return true;
}
vtkSmartPointer<vtkPolyData> TraceObject::GetVTKPolyData()
{
  
  hashp.clear();
  hashc.clear();
  //printf("Started creating vtkPolyData for rendering purposes ... ");
  vtkSmartPointer<vtkPolyData> poly_traces=vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkFloatArray> point_scalars=vtkSmartPointer<vtkFloatArray>::New();
  point_scalars->SetNumberOfComponents(1);
  vtkSmartPointer<vtkPoints> line_points=vtkSmartPointer<vtkPoints>::New();
  line_points->SetDataTypeToDouble();
  vtkSmartPointer<vtkCellArray> line_cells=vtkSmartPointer<vtkCellArray>::New();
  //printf("Starting CreatePolyDataRecursive\n");
  for(unsigned int counter=0; counter<trace_lines.size(); counter++)
  {
    //printf("Calling CreatePolyDataRecursive %dth time\n",counter+1);
    CreatePolyDataRecursive(trace_lines[counter],point_scalars,line_points,line_cells);
  }
  printf("Finished CreatePolyDataRecursive\n");
  poly_traces->SetPoints(line_points);
  poly_traces->SetLines(line_cells);

  poly_traces->GetPointData()->SetScalars(point_scalars);
  //printf("Done\n");
  return poly_traces;
}

bool TraceObject::ReadFromRPIXMLFile(char * filename)
{
  cout << "Started reading from " << filename << endl;
  TiXmlDocument doc(filename);
  doc.LoadFile();
  TiXmlHandle docHandle( &doc );
  TiXmlElement* lineElement =
    docHandle.FirstChild("Trace").FirstChild("TraceLine").Element();
  TiXmlElement* bitElement;
  float lineID;
  int lineType, lineParent, bitID;
  double bitX, bitY, bitZ;
  TraceLine * tline;
  vtksys::hash_map<unsigned int, unsigned long long int> hash_load;
  while(lineElement)
  {
    if(lineElement->QueryFloatAttribute("ID", &lineID) != TIXML_SUCCESS)
    {
      cerr << "ERROR: TraceLine has no ID" << endl;
      return false;
    }
    if(lineElement->QueryIntAttribute("Type", &lineType) != TIXML_SUCCESS)
    {
      lineType = 1;
    }
    if(lineElement->QueryIntAttribute("Parent", &lineParent) != TIXML_SUCCESS)
    {
      lineParent = -1;
    }
    if(hash_load.count(lineID)>0)
    {
      tline = reinterpret_cast<TraceLine*>(hash_load[lineID]);
    }
    else
    {
      tline = new TraceLine();
      hash_load[lineID] = reinterpret_cast<unsigned long long int>(tline);
    }
    tline->SetId(lineID);
    tline->SetType(lineType);
    //tline->setTraceColor(1.0/tline->GetType());
    //lauren added 11.6.2009
    switch( tline->GetType() )
    {
    case 1:
        tline->setTraceColor(.75);//cyan
        break;
    case 3:
    case 4:
    case 5:
        tline->setTraceColor(.90);//blue
        break;
    case 7:
    case 0:
    default:
        tline->setTraceColor(.25); //yellow
    }
    //end of added code
    


    if(lineParent != -1)
    {
      TraceLine *tparent;
      if(hash_load.count(lineParent)==0)
      {
        tparent = new TraceLine();
        hash_load[lineID] = reinterpret_cast<unsigned long long int>(tparent);
      }
      else
      {
        tparent = reinterpret_cast<TraceLine*>(hash_load[lineParent]);
      }
      tline->SetParent(tparent);
      tparent->GetBranchPointer()->push_back(tline);
    }
    else
    {
      trace_lines.push_back(tline);
    }
    bitElement = lineElement->FirstChildElement("TraceBit");
    if(!bitElement)
    {
      cerr << "Failed to initialize bitElement" << endl;
    }
    while(bitElement)
    {
      if(bitElement->QueryIntAttribute("ID", &bitID) != TIXML_SUCCESS)
      {
        cerr << "ERROR: TraceBit missing ID" << endl;
        return false;
      }
      if(bitElement->QueryDoubleAttribute("X", &bitX) != TIXML_SUCCESS)
      {
        cerr << "ERROR: TraceBit missing X value" << endl;
        return false;
      }
      if(bitElement->QueryDoubleAttribute("Y", &bitY) != TIXML_SUCCESS)
      {
        cerr << "ERROR: TraceBit missing Y value" << endl;
        return false;
      }
      if(bitElement->QueryDoubleAttribute("Z", &bitZ) != TIXML_SUCCESS)
      {
        cerr << "ERROR: TraceBit missing Z value" << endl;
        return false;
      }
      TraceBit tbit;
      tbit.x = bitX;
      tbit.y = bitY;
      tbit.z = bitZ;
      tbit.id = bitID;
      tline->AddTraceBit(tbit);
      bitElement = bitElement->NextSiblingElement();
    }
    lineElement = lineElement->NextSiblingElement();
  }
  return true;
}

void TraceObject::CollectIdsRecursive(std::vector<int> ids, TraceLine* tline)
{
  ids.push_back(tline->GetId());
  for(unsigned int counter = 0; counter < tline->GetBranchPointer()->size(); counter++)
  {
    this->CollectIdsRecursive(ids,(*tline->GetBranchPointer())[counter]);
  }
}
int TraceObject::getNewLineId()
{
  std::vector<int> ids;
  for(unsigned int counter=0; counter< trace_lines.size(); counter++)
  {
    this->CollectIdsRecursive(ids, trace_lines[counter]);
  }
  int newId = this->trace_lines.size();
  std::sort(ids.begin(),ids.end());
  for(unsigned int counter=0; counter < ids.size(); counter++)
  {
    if(newId == ids[counter])
    {
      newId++; // guarantees uniqueness because of sorted array.
    }
  }
  return newId;

  //bool uniqueId = false;
  //std::vector<TraceLine*>::iterator itr;
  //while(!uniqueId)
  //  {
  //  newId++;
  //  uniqueId = true;
  //  for(itr = this->trace_lines.begin(); itr != this->trace_lines.end(); itr++)
  //    {
  //    if((*itr)->GetId() == newId)
  //      {
  //      uniqueId = false;
  //      break;
  //      }
  //    //this assumes that a TraceLine will only have either 0 or 2 children
  //    if((*itr)->GetBranchPointer()->size() != 0)
  //      {
  //  for(int counter=0; counter < (*itr)->GetBranchPointer()->size(); counter++)
  //  {
  //    if( (*(*itr)->GetBranchPointer())[counter]->GetId() == newId)
  //    {
  //      uniqueId = false;
  //      break;
  //    }
  //  }
  //  if(uniqueId == false)
  //    break;
  //      }
  //    }
  //  }
  //return newId;
} 

/** For now splitTrace makes the following assumptions:
 *  1.  The newly created line should have the same type as the line it's being
 *      split from.
 *  2.  The newly created line is a "root line", ie it has no parent.
 *  Also note that we're automatically generating a new Id based on the number
 *  of trace lines that currently exist.  We test the uniqueness of the new ID
 *  before assignment.
 **/

void TraceObject::splitTrace(int selectedCellId)
{
  //get the TraceLine that contains the selected point
  TraceLine *selectedLine = 
    reinterpret_cast<TraceLine*>(this->hashc[selectedCellId]);

  if(selectedLine->GetSize() < 2)
  {
    cerr << "Cannot split a TraceLine that consists of fewer than two points, ";
    cerr << "Call (d)elete instead." << endl;
    return;
  }

  //initialize the new line that is being created by this split operation
  TraceLine *newLine = new TraceLine();
  newLine->SetType(selectedLine->GetType());
  int newId = this->getNewLineId();
  newLine->SetId(newId);

  //find the TraceBit that corresponds to selectedCellId
  std::list<TraceBit>::iterator bitItr =
    selectedLine->GetTraceBitIteratorBegin();
  std::vector<unsigned int>::iterator markerItr =
    selectedLine->GetMarkers()->begin(); 
  for(; markerItr != selectedLine->GetMarkers()->end() && bitItr != selectedLine->GetTraceBitIteratorEnd(); markerItr++)
    {
    if(*markerItr == (unsigned int)selectedCellId)
      {
      break;
      }
    bitItr++;
    }

  //some TraceLines have an equal number of cells and markers, whereas
  //other lines have one extra cell.  We have to treat these cases separately
  //to achieve consistent results.
  if((unsigned int)selectedLine->GetSize() > selectedLine->GetMarkers()->size())
    {
    bitItr++;
    }

  //special case: when splitting on the first node of a child line,
  //we need to remove selectedLine from the parent's vector of children.
  bool deleteSelectedLine = false;
  if( bitItr == selectedLine->GetTraceBitIteratorBegin() &&
      selectedLine->GetParent() != NULL )
    {
    //we want to do the actual removal after the splice, but we need to detect
    //this case before selectedLine->GetTraceBitIteratorBegin() changes
    deleteSelectedLine = true;
    }

  //move the TraceBits from selectedCellId on to the new TraceLine
  newLine->GetTraceBitsPointer()->splice(
    newLine->GetTraceBitIteratorBegin(),
    *selectedLine->GetTraceBitsPointer(),
    bitItr,
    selectedLine->GetTraceBitIteratorEnd());
  
  //if the selected line had any branches we have to reassign them to
  //the new line
  if(selectedLine->GetBranchPointer()->size() != 0)
    {
    *(newLine->GetBranchPointer()) = *(selectedLine->GetBranchPointer());
    std::vector<TraceLine*> * bp = newLine->GetBranchPointer();
    for(unsigned int counter=0; counter< bp->size(); counter++)
      {
      (*bp)[counter]->SetParent(newLine);
      }
    selectedLine->GetBranchPointer()->clear();
    }

  this->trace_lines.push_back(newLine);

  if(deleteSelectedLine)
    {
    std::vector<TraceLine*>::iterator tlitr = 
      find(selectedLine->GetParent()->GetBranchPointer()->begin(),
           selectedLine->GetParent()->GetBranchPointer()->end(),
           selectedLine);
    selectedLine->GetParent()->GetBranchPointer()->erase(tlitr);
    delete selectedLine;
    }
}

//This function assumes that the TraceLine does not have either parent or
//children. One end must be open 
void TraceObject::ReverseSegment(TraceLine *tline)
{
  if(tline->GetParent()==NULL)
  {
    if(tline->GetBranchPointer()->size()==0)
    {
      tline->GetTraceBitsPointer()->reverse();
      reverse(tline->GetMarkers()->begin(),tline->GetMarkers()->end());
      return;
    }
    else
    {
      tline->GetTraceBitsPointer()->reverse();
      reverse(tline->GetMarkers()->begin(),tline->GetMarkers()->end());
      // choose branch1 as parent arbitrarily
      TraceLine* temp1 =  tline->GetBranch1();
      TraceLine* temp2 = tline->GetBranch2();
      temp1->SetParent(NULL);
      ReverseSegment(temp1);
      tline->GetBranchPointer()->clear();
      tline->SetParent(temp1);
      temp2->SetParent(temp1);
      temp1->AddBranch(tline);
      temp1->AddBranch(tline);
    }
  }
  else if(tline->GetBranchPointer()->size()==0)
  {
    TraceLine * temp1 = tline->GetParent();
    TraceLine * temp2 = tline->GetParent()->GetBranch1();
    if(temp2 == tline)
    {
      temp2 = tline->GetParent()->GetBranch2();
    }
    tline->GetTraceBitsPointer()->reverse();
    reverse(tline->GetMarkers()->begin(),tline->GetMarkers()->end());
    temp1->GetBranchPointer()->clear();
    ReverseSegment(temp1);
    temp2->SetParent(tline);
    tline->SetParent(NULL);
    temp1->SetParent(tline);
    tline->AddBranch(temp1);
    tline->AddBranch(temp2);
  }
  else
  {
    //error!!
    printf("Error! the segment has to be open at one end before it can be reversed\n");
    return;
  }
}
//   t1, t2 -> two lines to be merged
//   pMarker - point of merging.

void TraceObject::FixPointMarkers(TraceLine* tline)
{
  TraceLine::TraceBitsType::iterator iter1 = tline->GetTraceBitIteratorBegin();
  TraceLine::TraceBitsType::iterator iter2 = tline->GetTraceBitIteratorEnd();
  while(iter1!=iter2)
  {
    hashp[iter1->marker]=reinterpret_cast<unsigned long long int>(tline);
    ++iter1;
  }
}
void TraceObject::mergeTraces(unsigned long long int eMarker, unsigned long long int sMarker)// point marker where we need to make the connection
{
  TraceLine * tmarker = reinterpret_cast<TraceLine*>(hashp[eMarker]);
  TraceLine * tother = reinterpret_cast<TraceLine*>(hashp[sMarker]);
  char elocation = -1;
  if(eMarker == tmarker->GetTraceBitsPointer()->front().marker)
  {
    elocation = 0;
  }
  else if(eMarker == tmarker->GetTraceBitsPointer()->back().marker)
  {
    elocation = 1;
  }

  char slocation = -1;
  if(sMarker == tother->GetTraceBitsPointer()->front().marker)
  {
    slocation = 0;
  }
  else if(sMarker == tother->GetTraceBitsPointer()->back().marker)
  {
    slocation = 1;
  }

  if(slocation ==0 && elocation ==1)
  {
    TraceLine::TraceBitsType::iterator iter = tother->GetTraceBitIteratorBegin();
    tmarker->GetTraceBitsPointer()->splice(tmarker->GetTraceBitIteratorEnd(),*(tother->GetTraceBitsPointer()));
    FixPointMarkers(tmarker);
    *(tmarker->GetBranchPointer())=*(tother->GetBranchPointer());
    tother->GetBranchPointer()->clear();
    for(unsigned int counter=0; counter< tmarker->GetBranchPointer()->size(); counter++)
    {
      (*tmarker->GetBranchPointer())[counter]->SetParent(tmarker);
    }
    RemoveTraceLine(tother);
    tmarker->setTraceColor(this->mergeLineColor); //tline->setTraceColor(1.0/tline->GetType());
  }
  else if (slocation ==1 && elocation == 0)
  {
    tother->GetTraceBitsPointer()->splice(tother->GetTraceBitIteratorEnd(),*(tmarker->GetTraceBitsPointer()));
    FixPointMarkers(tother);
    *(tother->GetBranchPointer())=*(tmarker->GetBranchPointer());
    tmarker->GetBranchPointer()->clear();
    for(unsigned int counter=0; counter< tother->GetBranchPointer()->size(); counter++)
    {
      (*tother->GetBranchPointer())[counter]->SetParent(tother);
    }
    RemoveTraceLine(tmarker);
    tother->setTraceColor(this->mergeLineColor);
  }
  else if (slocation == 0 && elocation ==0)
  {
    ReverseSegment(tother);
    tother->GetTraceBitsPointer()->splice(tother->GetTraceBitIteratorEnd(),*(tmarker->GetTraceBitsPointer()));
    FixPointMarkers(tother);
    *(tother->GetBranchPointer())=*(tmarker->GetBranchPointer());
    tmarker->GetBranchPointer()->clear();
    for(unsigned int counter=0; counter< tother->GetBranchPointer()->size(); counter++)
    {
      (*tother->GetBranchPointer())[counter]->SetParent(tother);
    }
    RemoveTraceLine(tmarker);
    tother->setTraceColor(this->mergeLineColor);//.7/tother->GetType()
  }
  else if (slocation == 1 && elocation ==1)
  {
    ReverseSegment(tother);
    tmarker->GetTraceBitsPointer()->splice(tmarker->GetTraceBitIteratorEnd(),*(tother->GetTraceBitsPointer()));
    FixPointMarkers(tmarker);
    *(tmarker->GetBranchPointer())=*(tother->GetBranchPointer());
    tother->GetBranchPointer()->clear();
    for(unsigned int counter=0; counter< tmarker->GetBranchPointer()->size(); counter++)
    {
      (*tmarker->GetBranchPointer())[counter]->SetParent(tmarker);
    }
    RemoveTraceLine(tother);
    tmarker->setTraceColor(this->mergeLineColor);
  }
  else if (slocation == -1 && elocation !=-1)
  {
std::cout<<" -1 !-1 "<<std::endl;
  }
  else if (slocation != -1 && elocation == -1)
  {
std::cout<<"!-1 -1"<<std::endl;
  }
  else
  {
    std::cout<<"else"<<std::endl;
  }


  //if(tother->GetParent()!=NULL)
  //{
  //  assert(tother->GetBranchPointer()->size()==0);//atleast one end must be open
  //  switch(plocation)
  //  {
  //  case 0:
  //    //Connect tother's back to tmarker's front

  //    break;
  //  case 1:
  //    //Connect tother's back to tmarker's back
  //    break;
  //  default:
  //    //Connect tother's back to tmarker's mid
  //    break;
  //  }
  //}
  //else if(tother->GetBranchPointer()->size()!=0)
  //{
  //  
  //  switch(plocation)
  //  {
  //  case 0:
  //    break;
  //  case 1:
  //    break;
  //  default:
  //    break;
  //  }
  //}
  //else
  //{
  //  TraceBit f = tother->GetTraceBitsPointer()->front();
  //  TraceBit b = tother->GetTraceBitsPointer()->back();
  //  TraceLine::TraceBitsType::iterator bit_iter = tmarker->GetTraceBitIteratorBegin();
  //  TraceLine::TraceBitsType::iterator bit_end = tmarker->GetTraceBitIteratorEnd();
  //  while(bit_iter!=bit_end)
  //  {
  //    if(bit_iter->marker==pMarker)
  //    {
  //      break;
  //    }
  //    ++bit_iter;
  //  }
  //  double dist1 = sqrt((f.x-bit_iter->x)*(f.x-bit_iter->x)+(f.y-bit_iter->y)*(f.y-bit_iter->y)+(f.z-bit_iter->z)*(f.z-bit_iter->z));
  //  double dist2 = sqrt((b.x-bit_iter->x)*(b.x-bit_iter->x)+(b.y-bit_iter->y)*(b.y-bit_iter->y)+(b.z-bit_iter->z)*(b.z-bit_iter->z));


  //  if(dist1 > dist2
  //}
  
}
void TraceObject::CollectBranchPointsRecursive(vtkSmartPointer<vtkPoints> p, vtkSmartPointer<vtkCellArray> cells,TraceLine *tline)
{
  if(tline->GetBranchPointer()->size()>0)
  {
    double loc[3];
    vtkIdType id;
    
    loc[0] = tline->GetTraceBitsPointer()->back().x;
    loc[1] = tline->GetTraceBitsPointer()->back().y;
    loc[2] = tline->GetTraceBitsPointer()->back().z;
    id = p->InsertNextPoint(loc);
    cells->InsertNextCell(1);
    cells->InsertCellPoint(id);
    for(unsigned int counter=0; counter< tline->GetBranchPointer()->size(); counter++)
    {
      this->CollectBranchPointsRecursive(p,cells,(*tline->GetBranchPointer())[counter]);
    }
  }
  
  return;
}
void TraceObject::CollectSegmentMidPointsRecursive(vtkSmartPointer<vtkPoints>p, vtkSmartPointer<vtkCellArray> cells, 
                    vtkSmartPointer<vtkFloatArray> da,TraceLine* tline)
{
    double loc[3];
    float dir[3];
    loc[0] = (tline->GetTraceBitsPointer()->back().x+tline->GetTraceBitsPointer()->front().x)/2.0;
    loc[1] = (tline->GetTraceBitsPointer()->back().y+tline->GetTraceBitsPointer()->front().y)/2.0;
    loc[2] = (tline->GetTraceBitsPointer()->back().z+tline->GetTraceBitsPointer()->front().z)/2.0;
    dir[0] = tline->GetTraceBitsPointer()->front().x-loc[0];
    dir[1] = tline->GetTraceBitsPointer()->front().y-loc[1];
    dir[2] = tline->GetTraceBitsPointer()->front().z-loc[2];

    if(tline->GetTraceBitsPointer()->size()==1)
    {
      if(tline->GetParent()!=NULL)
      {
        loc[0] = (tline->GetTraceBitsPointer()->back().x + tline->GetParent()->GetTraceBitsPointer()->back().x)/2.0;
        loc[1] = (tline->GetTraceBitsPointer()->back().y + tline->GetParent()->GetTraceBitsPointer()->back().y)/2.0;
        loc[2] = (tline->GetTraceBitsPointer()->back().z + tline->GetParent()->GetTraceBitsPointer()->back().z)/2.0;
        dir[0] = tline->GetParent()->GetTraceBitsPointer()->back().x-loc[0];
        dir[1] = tline->GetParent()->GetTraceBitsPointer()->back().y-loc[1];
        dir[2] = tline->GetParent()->GetTraceBitsPointer()->back().z-loc[2];
      }
    }
    vtkIdType id = p->InsertNextPoint(loc);
    cells->InsertNextCell(1);
    cells->InsertCellPoint(id);
    da->InsertNextTuple3(dir[0],dir[1],dir[2]);
    for(unsigned int counter=0; counter< tline->GetBranchPointer()->size(); counter++)
    {
      this->CollectSegmentMidPointsRecursive(p,cells,da,(*tline->GetBranchPointer())[counter]);
    }
  
}
vtkSmartPointer<vtkPolyData> TraceObject::generateBranchIllustrator()
{
  vtkSmartPointer<vtkSphereSource> s_src = vtkSmartPointer<vtkSphereSource>::New();
  s_src->SetRadius(1);

  VTK_CREATE(vtkArrowSource, arrow_src);

  VTK_CREATE(vtkGlyph3D, glyphs1);
  VTK_CREATE(vtkGlyph3D, glyphs2);
  VTK_CREATE(vtkPoints, p);
  VTK_CREATE(vtkCellArray, cells);
  

  printf("TraceLines size = %d\n",(int)trace_lines.size());
  for(unsigned int counter=0; counter< trace_lines.size(); counter++)
  {
    this->CollectBranchPointsRecursive(p,cells,trace_lines[counter]);
  }
  VTK_CREATE(vtkPoints, p1);
  VTK_CREATE(vtkFloatArray, da);
  da->SetNumberOfComponents(3);
  VTK_CREATE(vtkCellArray, cells1);
  for(unsigned int counter=0; counter< trace_lines.size(); counter++)
  {
    this->CollectSegmentMidPointsRecursive(p1,cells1,da,trace_lines[counter]);
  }

  VTK_CREATE(vtkPolyData, poly);
  poly->SetPoints(p);
  poly->SetVerts(cells);

  VTK_CREATE(vtkPolyData, poly1);
  poly1->SetPoints(p1);
  poly1->SetVerts(cells1);
  poly1->GetPointData()->SetVectors(da);

  glyphs2->SetInput(poly1);
  glyphs2->SetSource(arrow_src->GetOutput());
  glyphs2->SetVectorModeToUseVector();
  glyphs2->SetScaleFactor(5);

  glyphs1->SetInput(poly);
  glyphs1->SetSource(s_src->GetOutput());
  glyphs1->Update();
  glyphs2->Update();

  VTK_CREATE(vtkAppendPolyData, app_poly);
  app_poly->AddInput(glyphs1->GetOutput());
  app_poly->AddInput(glyphs2->GetOutput());
  app_poly->Update();
  
  return app_poly->GetOutput();
}

void TraceObject::RemoveTraceLine(TraceLine *tline)
{
  std::vector<TraceLine*>::iterator iter = trace_lines.begin();

  while(iter!=trace_lines.end())
  {
    //printf("TraceLine* = %p\n",*iter);
    if(tline == *iter)
    {
      trace_lines.erase(iter);
      return;
    }
    ++iter;
  }
//printf("Quitting RemoveTraceLine\n");
}

void TraceObject::FindMinLines(int smallSize)
{
  //std::cout<< "finding small lines\n";
  TraceLine *tline;
  std::vector<TraceLine*>::iterator iter = trace_lines.begin();
  while(iter!=trace_lines.end())
  {
    tline=*iter;
    if(smallSize >= tline->GetSize())
    {
      tline->setTraceColor(this->smallLineColor);
      SmallLines.push_back(tline);
    }
    ++iter;
  }
}

int TraceObject::createGapLists(std::vector<TraceLine*> traceList)
{ 
  unsigned int i,j, exist = 0, conflict = 0;  
  QProgressDialog progress("Searching for traces to merge",
                            "Cancel", 0, traceList.size() - 1);
  progress.setWindowModality(Qt::WindowModal);
  for (i=0;i<traceList.size()-1; i++)
    {
    progress.setValue(i);
    if(progress.wasCanceled())
      {
      return -1;
      }
    for (j=i+1; j<traceList.size(); j++)
      {
      TraceGap *newGap = new TraceGap;
      newGap->Trace1 = traceList[i];
      newGap->Trace2 = traceList[j];
      newGap->Trace1->EndPtDist(
          newGap->Trace2,newGap->endPT1, newGap->endPT2, 
          newGap->dist, newGap->maxdist, newGap->angle);
      newGap->length = newGap->Trace1->GetSize() + newGap->Trace2->GetSize() + newGap->dist;
      newGap->smoothness = newGap->length / newGap->maxdist;
      newGap->cost = newGap->angle*(newGap->dist/gapMax)*newGap->smoothness;
      if(!(newGap->dist >= newGap->Trace1->GetSize()*gapTol) &&
         !(newGap->dist >= newGap->Trace2->GetSize()*gapTol) &&
         !(newGap->dist >= gapMax*( 1+ gapTol)))
        { //myText+="added comparison\n";
        this->Gaps.push_back(newGap);
        } //end if
      }//end for j
    }// end for i
  if (this->Gaps.size() > 1)
    {   
    i = 0, j = 0;
    while (i < this->Gaps.size() -1)
      { //search for conflicts
      exist = 0;
      while ((exist == 0)&&(j<this->Gaps.size()-1))
        {
        j++;
        if (this->Gaps[i]->Trace1->GetId()==this->Gaps[j]->Trace1->GetId())
          {
          if (this->Gaps[i]->endPT1==this->Gaps[j]->endPT1)
            {
            exist = 1;
            }
          }
        else if(this->Gaps[i]->Trace1->GetId()==this->Gaps[j]->Trace2->GetId())
          {
          if (this->Gaps[i]->endPT1==this->Gaps[j]->endPT2)
            {
            exist = 1;
            }
          }
        else if (this->Gaps[i]->Trace2->GetId() == this->Gaps[j]->Trace1->GetId())
          {
          if (this->Gaps[i]->endPT2==this->Gaps[j]->endPT1)
            {
            exist = 1;
            }
          }
        else if(this->Gaps[i]->Trace2->GetId() == this->Gaps[j]->Trace2->GetId())
          {
          if (this->Gaps[i]->endPT2==this->Gaps[j]->endPT2)
            {
            exist = 1;
            }
          }
        }   //end while exist = 0
      if (exist == 1)
        {
        ++conflict;
        if (this->Gaps[i]->cost<this->Gaps[j]->cost)
          {
          this->Gaps.erase(this->Gaps.begin()+j);
          }
        else
          {
          this->Gaps.erase(this->Gaps.begin()+i);
          }
        j = i;
        }//end if exist
      else
        {
        i++;
        j=i;
        }//end else exist
      }// end of search for conflicts
    }
  return conflict;
}
