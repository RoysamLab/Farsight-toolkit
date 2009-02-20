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
  point[0] = iter->x;point[1]=iter->y;point[2]=iter->z;
  return_id = line_points->InsertNextPoint(point);
  hashp[return_id]=(unsigned long long int)tline;
  iter->marker = return_id;
  point_scalars->InsertNextTuple1(1.0/tline->GetType());

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

    point_scalars->InsertNextTuple1(1.0/tline->GetType());
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

struct equllint
{
  bool operator()(const unsigned long long int i1, unsigned long long int i2) const
    {
    return i1 == i2;
    }
};

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
};


bool TraceObject::WriteToSWCFile(char *filename)
{
  FILE * fp = fopen(filename,"w");
  hash_map<const unsigned long long int, int, hashulli> hash_dump;
  //hash_map<const unsigned long long int, int, hashulli, equlli > hash_dump;
  //hash_map<const unsigned long long int, int> hash_dump;
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
    hash_dump[reinterpret_cast<const unsigned long long int>(t)]=cur_id+t->GetTraceBitsPointer()->size()-1;
    if(t->GetParent()==NULL)
    {
      fprintf(fp,"%d %d %0.2lf %0.2lf %0.2lf %0.2lf %d\n",cur_id++,t->GetType(),iter->x,iter->y,iter->z,iter->r,-1);
    }
    else
    {
      fprintf(fp,"%d %d %0.2lf %0.2lf %0.2lf %0.2lf %d\n",cur_id++,t->GetType(),iter->x,iter->y,iter->z,iter->r,hash_dump[reinterpret_cast<const unsigned long long int>(t->GetParent())]);
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
    printf("Calling CreatePolyDataRecursive %dth time\n",counter+1);
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
  if(!lineElement)
    {
    cerr << "Problem chief: failed to initialize lineElement" << endl;
    }
  while(lineElement)
    {
    if(lineElement->QueryFloatAttribute("ID", &lineID) != TIXML_SUCCESS)
      {
      cerr << "ERROR: TraceLine has no ID" << endl;
      return false;
      }
    cout << "Found line ID #" << lineID << endl;
    if(lineElement->QueryIntAttribute("Type", &lineType) != TIXML_SUCCESS)
      {
      lineType = 1;
      }
    if(lineElement->QueryIntAttribute("Parent", &lineParent) != TIXML_SUCCESS)
      {
      lineParent = -1;
      }
    tline = new TraceLine();
    tline->SetId(lineID);
    tline->SetType(lineType);
    hash_load[lineID] = reinterpret_cast<unsigned long long int>(tline);
    if(lineParent != -1)
      {
      TraceLine * tparent = reinterpret_cast<TraceLine*>(hash_load[lineParent]);
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
      cout << "Found bit ID #" << bitID << endl;
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
}
