#include "Trace.h"
#pragma warning(disable:4996)

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
  hash_map<unsigned int, unsigned long long int> hash_load;
  int tc = 0;
  unsigned char child_count[1000000];
      int id, type,parent;
    double x,y,z,r;
  while(!feof(fp))
  {
    fgets(buff,1024,fp);
    int pc = 0;
    while(buff[pc]==' '&&(pc<1023))
      pc++;
    if(buff[pc]=='#') // ignoring comment lines for now. We have to read scale from it! TODO
      continue;

    //sscanf(buff,"%d %d %lf %lf %lf %lf %d",&id,&type,&x,&y,&z,&r,&parent);
    sscanf(buff,"%d %*lf %*lf %*lf %*lf %*lf %d",&id,&parent);
    tc++;
    child_count[id]=0;
    if(parent!=-1)
      child_count[parent]++;
  }
  //printf("I read %d lines\n",tc);

  fclose(fp);
  fp = fopen(filename,"r");
  int tcc =0;
  while(!feof(fp))
  {
    tcc++;
    //printf("Done %d\n",tcc);
    fgets(buff,1024,fp);
    int pc = 0;
    while(buff[pc]==' '&&(pc<1023))
      pc++;
    if(buff[pc]=='#') // ignoring comment lines for now. We have to read scale from it! TODO
      continue;
    sscanf(buff,"%d %d %lf %lf %lf %lf %d",&id,&type,&x,&y,&z,&r,&parent);
    TraceBit tbit;
    tbit.x=x;tbit.y=y;tbit.z=z;tbit.id=id;tbit.r =r;
    if(parent==-1)
    {
      TraceLine *line = new TraceLine();  
      line->SetType(type);
      line->SetId(type);
      line->AddTraceBit(tbit);
      hash_load[id]=(unsigned long long int)line;
      trace_lines.push_back(line);
    }
    else if (child_count[parent]==1)
    {
        
        TraceLine* line = reinterpret_cast<TraceLine*>(hash_load[parent]);
        hash_load[id] = (unsigned long long int)line;
        TraceLine::TraceBitsType::iterator iter = line->GetTraceBitIteratorEnd();
        iter--;
        TraceLine::TraceBitsType::iterator iterbegin = line->GetTraceBitIteratorBegin();
        do
        {
          if(iter->id==parent)
          {
            ++iter;
            line->GetTraceBitsPointer()->insert(iter,tbit);
          //  printf("Successfully inserted\n");
            break;
          }
          --iter;
        }while(iter!=iterbegin);
    }
    else 
    {
      TraceLine * line = new TraceLine();
      line->SetType(type);
      line->SetId(type);
      line->AddTraceBit(tbit);
      hash_load[id]=(unsigned long long int)line;
      TraceLine * old_line = reinterpret_cast<TraceLine*>(hash_load[parent]);
      old_line->AddBranch(line);
      line->SetParent(old_line);
    }
  }
  printf("Finished loading\n");
  //Print(std::cout);
  fclose(fp);
  return true;
}

void TraceObject::CreatePolyDataRecursive(TraceLine* tline, vtkSmartPointer<vtkFloatArray> point_scalars, vtkSmartPointer<vtkPoints> line_points,vtkSmartPointer<vtkCellArray> line_cells)
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
  point_scalars->InsertNextTuple1(1-1.0/tline->GetType());

  /* To add a line between parent line's last point and the first point in the current line */
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
  /* Rest of the lines for the current tline */
  iter++;
  int pc = 0;
  while(iter!=tline->GetTraceBitIteratorEnd())
  {
    //printf("in loop %d\n",++pc);
    old_id = return_id;
    point[0] = iter->x;point[1]=iter->y;point[2]=iter->z;
    return_id = line_points->InsertNextPoint(point);
    hashp[return_id]=(unsigned long long int)tline;
    iter->marker = return_id;

    point_scalars->InsertNextTuple1(1-1.0/tline->GetType());
    cell_id = line_cells->InsertNextCell(2);
    cell_id_array->push_back(cell_id);
    hashc[cell_id]=reinterpret_cast<unsigned long long int>(tline);
    line_cells->InsertCellPoint(old_id);
    line_cells->InsertCellPoint(return_id);
    ++iter;
  }
  /* Recursive calls to the branches if they exist */
  for(int counter=0; counter< tline->GetBranchPointer()->size(); counter++)
  {
    //printf("I should be having children too! what am I doing here?\n");
    CreatePolyDataRecursive((*tline->GetBranchPointer())[counter],point_scalars,line_points,line_cells);
  }
}

void CollectTraceBitsRecursive(std::vector<TraceBit> &vec,TraceLine *l)
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
  for(int counter=0; counter< bp->size(); counter++)
  {
    CollectTraceBitsRecursive(vec,(*bp)[counter]);
  }
}

std::vector<TraceBit> TraceObject::CollectTraceBits()
{
  std::vector<TraceBit> vec;
  for(int counter=0; counter<trace_lines.size(); counter++)
    CollectTraceBitsRecursive(vec,trace_lines[counter]);
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

bool TraceObject::WriteToSWCFile(char *filename)
{
  FILE * fp = fopen(filename,"w");
  hash_map<const unsigned long long int, int, hashulli> hash_dump;
  if(fp == NULL)
  {
    printf("Couldn't open %s for writing\n",filename);
    return false;
  }
  int cur_id = 1;
  std::queue<TraceLine*> q;
  for(int counter=0; counter<trace_lines.size(); counter++)
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
    for(int counter=0; counter<t->GetBranchPointer()->size(); counter++)
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
  printf("Started creating vtkPolyData for rendering purposes ... ");
  vtkSmartPointer<vtkPolyData> poly_traces=vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkFloatArray> point_scalars=vtkSmartPointer<vtkFloatArray>::New();
  point_scalars->SetNumberOfComponents(1);
  vtkSmartPointer<vtkPoints> line_points=vtkSmartPointer<vtkPoints>::New();
  line_points->SetDataTypeToDouble();
  vtkSmartPointer<vtkCellArray> line_cells=vtkSmartPointer<vtkCellArray>::New();
  printf("Starting CreatePolyDataRecursive\n");
  for(int counter=0; counter<trace_lines.size(); counter++)
  {
    /*printf("Calling CreatePolyDataRecursive %dth time\n",counter+1);*/
    CreatePolyDataRecursive(trace_lines[counter],point_scalars,line_points,line_cells);
  }
  printf("Finished CreatePolyDataRecursive\n");
  poly_traces->SetPoints(line_points);
  poly_traces->SetLines(line_cells);
  
  poly_traces->GetPointData()->SetScalars(point_scalars);
  printf("Done\n");
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
  hash_map<unsigned int, unsigned long long int> hash_load;
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

void CollectIdsRecursive(std::vector<int> ids, TraceLine* tline)
{
	ids.push_back(tline->GetId());
	for(int counter = 0; counter < tline->GetBranchPointer()->size(); counter++)
	{
		CollectIdsRecursive(ids,(*tline->GetBranchPointer())[counter]);
	}
}
int TraceObject::getNewLineId()
{
  std::vector<int> ids;
  for(int counter=0; counter< trace_lines.size(); counter++)
  {
	  CollectIdsRecursive(ids, trace_lines[counter]);
  }
  int newId = this->trace_lines.size();
  std::sort(ids.begin(),ids.end());
  for(int counter=0; counter < ids.size(); counter++)
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
		//	for(int counter=0; counter < (*itr)->GetBranchPointer()->size(); counter++)
		//	{
		//		if( (*(*itr)->GetBranchPointer())[counter]->GetId() == newId)
		//		{
		//			uniqueId = false;
		//			break;
		//		}
		//	}
		//	if(uniqueId == false)
		//		break;
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
    bitItr++;
    if(*markerItr == selectedCellId)
      {
      break;
      }
    }

  //move the TraceBits from selectedCellId on to the new TraceLine
  newLine->GetTraceBitsPointer()->splice(
    newLine->GetTraceBitIteratorBegin(),
    *selectedLine->GetTraceBitsPointer(),
    bitItr,
    selectedLine->GetTraceBitIteratorEnd());
  
  //if the selected line had any branches, we have to reassign them to
  //the new line
  if(selectedLine->GetBranchPointer()->size() != 0)
    {
    *(newLine->GetBranchPointer()) = *(selectedLine->GetBranchPointer());
	
	std::vector<TraceLine*> * bp = newLine->GetBranchPointer();
	for(int counter=0; counter< bp->size(); counter++)
	{
		(*bp)[counter]->SetParent(newLine);
	}
    selectedLine->GetBranchPointer()->clear();
   // newLine->GetBranch1()->SetParent(newLine);
   // newLine->GetBranch2()->SetParent(newLine);
    }

  this->trace_lines.push_back(newLine);
}

/* This function assumes that the TraceLine does not have either parent or children. One end must be open */
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
/*
t1, t2 -> two lines to be merged
pMarker - point of merging.



*/
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
		tmarker->GetTraceBitsPointer()->splice(tmarker->GetTraceBitIteratorEnd(),*(tother->GetTraceBitsPointer()));
		*(tmarker->GetBranchPointer())=*(tother->GetBranchPointer());
	}
	else if (slocation ==1 && elocation == 0)
	{
		tother->GetTraceBitsPointer()->splice(tother->GetTraceBitIteratorEnd(),*(tmarker->GetTraceBitsPointer()));
		*(tother->GetBranchPointer())=*(tmarker->GetBranchPointer());
	}
	else if (slocation == 0 && elocation ==0)
	{
		ReverseSegment(tother);
		tother->GetTraceBitsPointer()->splice(tother->GetTraceBitIteratorEnd(),*(tmarker->GetTraceBitsPointer()));
		*(tother->GetBranchPointer())=*(tmarker->GetBranchPointer());
	}
	else if (slocation == 1 && elocation ==1)
	{
		ReverseSegment(tother);
		tmarker->GetTraceBitsPointer()->splice(tmarker->GetTraceBitIteratorEnd(),*(tother->GetTraceBitsPointer()));
		*(tmarker->GetBranchPointer())=*(tother->GetBranchPointer());
	}
	else if (slocation == -1 && elocation !=-1)
	{
		
	}
	else if (slocation != -1 && elocation == -1)
	{

	}
	else
	{

	}
	
	
	//if(tother->GetParent()!=NULL)
	//{
	//	assert(tother->GetBranchPointer()->size()==0);//atleast one end must be open
	//	switch(plocation)
	//	{
	//	case 0:
	//		//Connect tother's back to tmarker's front

	//		break;
	//	case 1:
	//		//Connect tother's back to tmarker's back
	//		break;
	//	default:
	//		//Connect tother's back to tmarker's mid
	//		break;
	//	}
	//}
	//else if(tother->GetBranchPointer()->size()!=0)
	//{
	//	
	//	switch(plocation)
	//	{
	//	case 0:
	//		break;
	//	case 1:
	//		break;
	//	default:
	//		break;
	//	}
	//}
	//else
	//{
	//	TraceBit f = tother->GetTraceBitsPointer()->front();
	//	TraceBit b = tother->GetTraceBitsPointer()->back();
	//	TraceLine::TraceBitsType::iterator bit_iter = tmarker->GetTraceBitIteratorBegin();
	//	TraceLine::TraceBitsType::iterator bit_end = tmarker->GetTraceBitIteratorEnd();
	//	while(bit_iter!=bit_end)
	//	{
	//		if(bit_iter->marker==pMarker)
	//		{
	//			break;
	//		}
	//		++bit_iter;
	//	}
	//	double dist1 = sqrt((f.x-bit_iter->x)*(f.x-bit_iter->x)+(f.y-bit_iter->y)*(f.y-bit_iter->y)+(f.z-bit_iter->z)*(f.z-bit_iter->z));
	//	double dist2 = sqrt((b.x-bit_iter->x)*(b.x-bit_iter->x)+(b.y-bit_iter->y)*(b.y-bit_iter->y)+(b.z-bit_iter->z)*(b.z-bit_iter->z));


	//	if(dist1 > dist2
	//}
}
void CollectBranchPointsRecursive(vtkSmartPointer<vtkPoints> p, vtkSmartPointer<vtkCellArray> cells, vtkSmartPointer<vtkFloatArray> da, TraceLine *tline)
{
	if(tline->GetBranchPointer()->size()>0)
	{
		double loc[3];
		loc[0] = tline->GetTraceBitsPointer()->back().x;
		loc[1] = tline->GetTraceBitsPointer()->back().y;
		loc[2] = tline->GetTraceBitsPointer()->back().z;
		float dir[3];
		dir[0] = tline->GetTraceBitsPointer()->front().x-loc[0];
		dir[1] = tline->GetTraceBitsPointer()->front().y-loc[1];
		dir[2] = tline->GetTraceBitsPointer()->front().z-loc[2];
		vtkIdType id = p->InsertNextPoint(loc);
		cells->InsertNextCell(1);
		cells->InsertCellPoint(id);
		if(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]<1e-6)
		{
			printf("Error! Error! Error!\n");
			if(tline->GetParent()!=NULL)
			{
				dir[0] = -tline->GetTraceBitsPointer()->front().x+tline->GetParent()->GetTraceBitsPointer()->back().x;
				dir[1] = -tline->GetTraceBitsPointer()->front().y+tline->GetParent()->GetTraceBitsPointer()->back().y;
				dir[2] = -tline->GetTraceBitsPointer()->front().z+tline->GetParent()->GetTraceBitsPointer()->back().z;
			}
		}
		da->InsertNextTuple3(dir[0],dir[1],dir[2]);

		for(int counter=0; counter< tline->GetBranchPointer()->size(); counter++)
		{
			CollectBranchPointsRecursive(p,cells,da,(*tline->GetBranchPointer())[counter]);
		}
	}
	return;
}
vtkSmartPointer<vtkPolyData> TraceObject::generateBranchIllustrator()
{
	vtkSmartPointer<vtkSphereSource> s_src = vtkSmartPointer<vtkSphereSource>::New();
	s_src->SetRadius(1.0);

	VTK_CREATE(vtkArrowSource, arrow_src);

	VTK_CREATE(vtkGlyph3D, glyphs1);
	VTK_CREATE(vtkGlyph3D, glyphs2);
	VTK_CREATE(vtkPoints, p);
	VTK_CREATE(vtkFloatArray, da);
	VTK_CREATE(vtkCellArray, cells);
	da->SetNumberOfComponents(3);

	printf("TraceLines size = %d\n",trace_lines.size());
	for(int counter=0; counter< trace_lines.size(); counter++)
	{
		CollectBranchPointsRecursive(p,cells,da,trace_lines[counter]);
	}
	VTK_CREATE(vtkPolyData, poly);
	poly->SetPoints(p);
	poly->SetVerts(cells);
	poly->GetPointData()->SetVectors(da);

	glyphs2->SetInput(poly);
	glyphs2->SetSource(arrow_src->GetOutput());
	glyphs2->SetVectorModeToUseVector();
	glyphs2->SetScaleFactor(5);
	glyphs1->SetInput(poly);
	glyphs1->SetSource(s_src->GetOutput());
	glyphs1->Update();
	glyphs2->Update();
	return glyphs2->GetOutput();
}
void TraceLine::Getstats()
{
  double XF, XB, YF, YB, ZF, ZB;
  XF = m_trace_bits.front().x;
  XB = m_trace_bits.back().x;
  YF= m_trace_bits.front().y;
  YB= m_trace_bits.back().y;
  ZF= m_trace_bits.front().z;
  ZB= m_trace_bits.back().z;
  printf("Trace # %d \t Trace Size: \t %d ", m_id, m_trace_bits.size());
  printf("First bit x: %4.2f y: %4.2f z: %4.2f \t", XF, YF, ZF); 
  printf("Endt bit x: %4.2f y: %4.2f z: %4.2f \n", XB, YB, ZB); 

}
void TraceLine::EndPtDist(TraceLine *Trace2, int &dir1, int &dir2, double &dist)	
{
	//compute eucliedan distance 
	double XF, XB, YF, YB, ZF, ZB, XF2,XB2, YF2, YB2, ZF2, ZB2, distances[4], min ;

	XF = m_trace_bits.front().x;	YF= m_trace_bits.front().y;		ZF= m_trace_bits.front().z;
	XB = m_trace_bits.back().x;		YB= m_trace_bits.back().y;		ZB= m_trace_bits.back().z;
	XF2=Trace2->m_trace_bits.front().x;		YF2=Trace2->m_trace_bits.front().y;		ZF2=Trace2->m_trace_bits.front().z;
	XB2=Trace2->m_trace_bits.back().x;		YB2=Trace2->m_trace_bits.back().y;		ZB2=Trace2->m_trace_bits.back().z;
//compute the endpt distances
	distances[0]=sqrt(pow((XF-XF2),2)+pow((YF-YF2),2)+pow((ZF-ZF2),2));//0 F-F
	distances[1]=sqrt(pow((XF-XB2),2)+pow((YF-YB2),2)+pow((ZF-ZB2),2));//1 F-B
	
	distances[2]=sqrt(pow((XB-XF2),2)+pow((YB-YF2),2)+pow((ZB-ZF2),2));//2 B-F
	distances[3]=sqrt(pow((XB-XB2),2)+pow((YB-YB2),2)+pow((ZB-ZB2),2));//3 B-B
//determine minimum spacing
	min = distances[1];
	int i, mark=0;
	for (i = 2; i<4; i++)
	{
		if (min > distances[i])
		{
			min= distances[i];
			mark=i;
		}
	}
// from min determine orientation and return distance
	if (mark ==0)
	{
		dist = distances[0];
		dir1= m_trace_bits.front().id;
		dir2= Trace2->m_trace_bits.front().id;
	}
	else if (mark ==1)
	{
		dist = distances[1];
		dir1= m_trace_bits.front().id;
		dir2= Trace2->m_trace_bits.back().id;
	}
	else if (mark ==2)
	{
		dist = distances[2];
		dir1= m_trace_bits.back().id;
		dir2= Trace2->m_trace_bits.front().id;
	}
	else
	{
		dist = distances[3];
		dir1= m_trace_bits.back().id;
		dir2= Trace2->m_trace_bits.back().id;
	}
}

