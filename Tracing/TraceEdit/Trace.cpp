#include "Trace.h"
//#pragma warning(disable:4996)

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
	
	int tc = 0;
	const int MAXPOINTSINSWC= 10000;
	unsigned char child_count[MAXPOINTSINSWC];//MAX POINTS in the swc file Fix: TODO : probably make it a hash_map which is slower
	std::set<int> criticals; // store all points who have parent = -1 or child_count[parent] > 1

  vtksys::hash_map<unsigned int,int> hash_type; // smaller hash functions only for the critical points.. saves us memory and time
  vtksys::hash_map<unsigned int,int> hash_parent;
	vtksys::hash_map<unsigned int,unsigned long long int> hash_load;

	//memset(child_count,0,sizeof(unsigned char)*100000);
	
	for(int counter=0; counter < 10000; counter++)
		child_count[counter]=0;
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

		sscanf(buff,"%d %d %*f %*f %*f %*f %d",&id,&type,&parent);
		//sscanf(buff,"%d %d %*lf %*lf %*lf %*lf %d",&id,&type,&parent);
		tc++;
		//printf("%d\n",id);
		if(id>max_id)//find max id
		{
			max_id=id;
		}

		if(parent != -1)
		{
			child_count[parent]++;
		}
	}
	fclose(fp);
	//printf("I read %d lines\n",tc);
	unsigned int child_id[10000];
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
			hash_type[id] = type;
			hash_parent[id] = -1;
			//printf("hash_parent[%d] = %d\n",id,hash_parent[id]);
		}
		else
		{
			if(child_count[parent]>1)
			{
				criticals.insert(id);
				hash_type[id] = type;
				hash_parent[id] = parent;
				//printf("hash_parent[%d] = %d\n", id, hash_parent[id]);
			}
		}
	}
	
	std::set<int>::iterator iter = criticals.begin();
	int global_id_number = 1;
	while(iter != criticals.end())
	{
		TraceLine * ttemp = new TraceLine();
		ttemp->SetId(global_id_number++);
		ttemp->SetType(hash_type[*iter]);
		ttemp->setTraceColor(1.0/ttemp->GetType());
		//tline->setTraceColor(1.0/tline->GetType());
		ttemp->AddTraceBit(data[*iter]);
		int id_counter = *iter;
		while(child_count[id_counter]==1)
		{
			id_counter = child_id[id_counter];
			ttemp->AddTraceBit(data[id_counter]);
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
			//	printf("here is the error\n");
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
	fclose(fp);
//	delete [] child_id;
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
	point_scalars->InsertNextTuple1(.5-tline->getTraceColor());

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
	while(iter!=tline->GetTraceBitIteratorEnd())
	{
		//printf("in loop %d\n",++pc);
		old_id = return_id;
		point[0] = iter->x;point[1]=iter->y;point[2]=iter->z;
		return_id = line_points->InsertNextPoint(point);
		hashp[return_id]=(unsigned long long int)tline;
		iter->marker = return_id;

		point_scalars->InsertNextTuple1(tline->getTraceColor()-.25);
		cell_id = line_cells->InsertNextCell(2);
		cell_id_array->push_back(cell_id);
		hashc[cell_id]=reinterpret_cast<unsigned long long int>(tline);
		line_cells->InsertCellPoint(old_id);
		line_cells->InsertCellPoint(return_id);
		++iter;
	}
	/* Recursive calls to the branches if they exist */
	for(unsigned int counter=0; counter<tline->GetBranchPointer()->size(); counter++)
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
	for(unsigned int counter=0; counter< bp->size(); counter++)
	{
		CollectTraceBitsRecursive(vec,(*bp)[counter]);
	}
}

std::vector<TraceBit> TraceObject::CollectTraceBits()
{
	std::vector<TraceBit> vec;
	for(unsigned int counter=0; counter<trace_lines.size(); counter++)
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
		/*printf("Calling CreatePolyDataRecursive %dth time\n",counter+1);*/
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
		tline->setTraceColor(1.0/tline->GetType());
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
	for(unsigned int counter = 0; counter < tline->GetBranchPointer()->size(); counter++)
	{
		CollectIdsRecursive(ids,(*tline->GetBranchPointer())[counter]);
	}
}
int TraceObject::getNewLineId()
{
	std::vector<int> ids;
	for(unsigned int counter=0; counter< trace_lines.size(); counter++)
	{
		CollectIdsRecursive(ids, trace_lines[counter]);
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
		tmarker->setTraceColor(.7/tmarker->GetType()); //tline->setTraceColor(1.0/tline->GetType());
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
		tother->setTraceColor(.7/tother->GetType());
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
		tother->setTraceColor(.7/tother->GetType());
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
		tmarker->setTraceColor(.7/tmarker->GetType());
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
void CollectBranchPointsRecursive(vtkSmartPointer<vtkPoints> p, vtkSmartPointer<vtkCellArray> cells,TraceLine *tline)
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
			CollectBranchPointsRecursive(p,cells,(*tline->GetBranchPointer())[counter]);
		}
	}
	
	return;
}
void CollectSegmentMidPointsRecursive(vtkSmartPointer<vtkPoints>p, vtkSmartPointer<vtkCellArray> cells, vtkSmartPointer<vtkFloatArray> da,TraceLine* tline)
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
			CollectSegmentMidPointsRecursive(p,cells,da,(*tline->GetBranchPointer())[counter]);
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
		CollectBranchPointsRecursive(p,cells,trace_lines[counter]);
	}
	VTK_CREATE(vtkPoints, p1);
	VTK_CREATE(vtkFloatArray, da);
	da->SetNumberOfComponents(3);
	VTK_CREATE(vtkCellArray, cells1);
	for(unsigned int counter=0; counter< trace_lines.size(); counter++)
	{
		CollectSegmentMidPointsRecursive(p1,cells1,da,trace_lines[counter]);
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
void TraceLine::Getstats()
{
	double XF, XB, YF, YB, ZF, ZB;
	XF = m_trace_bits.front().x;
	XB = m_trace_bits.back().x;
	YF= m_trace_bits.front().y;
	YB= m_trace_bits.back().y;
	ZF= m_trace_bits.front().z;
	ZB= m_trace_bits.back().z;
	printf("Trace # %d \t Trace Size: \t %d ", m_id, (int)m_trace_bits.size());
	printf("First bit x: %4.2f y: %4.2f z: %4.2f \t", XF, YF, ZF); 
	printf("Endt bit x: %4.2f y: %4.2f z: %4.2f \n", XB, YB, ZB); 

}
void TraceLine::EndPtDist(TraceLine *Trace2, int &dir1, int &dir2, double &dist, double &maxdist, double &angle)	
{
	//compute eucliedan distance 
	double XF, XB, YF, YB, ZF, ZB, XF2,XB2, YF2, YB2, ZF2, ZB2, distances[4];
	double delX1, delX2, delY1, delY2, delZ1, delZ2, min, norm1, norm2;
//trace 1
	XF = m_trace_bits.front().x;	YF= m_trace_bits.front().y;		ZF= m_trace_bits.front().z;
	XB = m_trace_bits.back().x;		YB= m_trace_bits.back().y;		ZB= m_trace_bits.back().z;
//trace 2 
	XF2=Trace2->m_trace_bits.front().x;		YF2=Trace2->m_trace_bits.front().y;		ZF2=Trace2->m_trace_bits.front().z;
	XB2=Trace2->m_trace_bits.back().x;		YB2=Trace2->m_trace_bits.back().y;		ZB2=Trace2->m_trace_bits.back().z;
//delta x,y,z 
	delX1=XF-XB;	delX2=XF2-XB2;
	delY1=YF-YB;	delY2=YF2-YB2;
	delZ1=ZF-ZB;	delZ2=ZF2-ZB2;

	//compute the endpt distances
	distances[0]=sqrt(pow((XF-XF2),2)+pow((YF-YF2),2)+pow((ZF-ZF2),2));//0 F-F
	//std::cout<<distances[0]<<std::endl;
	distances[1]=sqrt(pow((XF-XB2),2)+pow((YF-YB2),2)+pow((ZF-ZB2),2));//1 F-B
	//std::cout<<distances[1]<<std::endl;
	distances[2]=sqrt(pow((XB-XF2),2)+pow((YB-YF2),2)+pow((ZB-ZF2),2));//2 B-F
	//std::cout<<distances[2]<<std::endl;
	distances[3]=sqrt(pow((XB-XB2),2)+pow((YB-YB2),2)+pow((ZB-ZB2),2));//3 B-B

	norm1=sqrt(pow((delX1),2)+ pow((delY1),2)+ pow((delZ1),2));
	norm2=sqrt(pow((delX2),2)+ pow((delY2),2)+ pow((delZ2),2));
	angle = acos(((delX1 * delX2) + (delY1 *delY2) + (delZ1 * delZ2))/(norm2 * norm1));
	//std::cout<<distances[3]<<std::endl;
	//determine minimum spacing
	min = distances[0];
	int i, mark=0;
	for (i = 1; i<4; i++)
	{
		//printf("dist of %d = %d \n", i, distances[i]);
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
		dir1= m_trace_bits.front().marker;
		dir2= Trace2->m_trace_bits.front().marker;
		maxdist= distances[3];
	}
	else if (mark ==1)
	{
		dist = distances[1];
		dir1= m_trace_bits.front().marker;
		dir2= Trace2->m_trace_bits.back().marker;
		maxdist= distances[2];
	}
	else if (mark ==2)
	{
		dist = distances[2];
		dir1= m_trace_bits.back().marker;
		dir2= Trace2->m_trace_bits.front().marker;
		maxdist= distances[1];
	}
	else
	{
		dist = distances[3];
		dir1= m_trace_bits.back().marker;
		dir2= Trace2->m_trace_bits.back().marker;
		maxdist= distances[0];
	}
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
			tline->setTraceColor(.3);
			SmallLines.push_back(tline);
		}
		++iter;
	}
}