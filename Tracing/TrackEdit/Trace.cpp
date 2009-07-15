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

#include "trace.h"
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

		tbit.x = x;	tbit.y = 2.97/2.79*y;	tbit.z = t;	tbit.id = z;line_count++;		
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
	stdext::hash_map<unsigned int, unsigned long long int> hash_load;
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
					//	printf("Successfully inserted\n");
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
	point_scalars->InsertNextTuple1(iter->id/40.0);

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

		point_scalars->InsertNextTuple1(iter->id/40.0);
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

bool TraceObject::WriteToSWCFile(char *filename)
{
	FILE * fp = fopen(filename,"w");
	stdext::hash_map<unsigned long long int,int> hash_dump;
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


// Handy functions to read Props from the xml
// Hope this is the right way to read and delete memory
//int xmlReadInt(xmlNodePtr p,char *s,int def)
//{
//	xmlChar * charp = xmlGetProp(p, BAD_CAST s);
//	int num;
//	if(charp!=NULL)
//	{
//		num = atof((char*)charp);
//		delete [] charp;
//	}
//	else
//	{
//		num = def; //	THIS SHOULD NEVER HAPPEN 
//	}
//	return num;
//}
//
//float xmlReadFloat(xmlNodePtr p,char *s,float def)
//{
//	xmlChar * charp = xmlGetProp(p, BAD_CAST s);
//	float num;
//	if(charp!=NULL)
//	{
//		num = atof((char*)charp);
//		delete [] charp;
//	}
//	else
//	{
//		num = def; //	THIS SHOULD NEVER HAPPEN 
//	}
//	return num;
//}
//
//bool TraceObject::ReadFromRPIXMLFile(char * filename)
//{
//	printf("Started reading from %s ...",filename);
//
//	stdext::hash_map<unsigned int, unsigned long long int> hash_load;
//	/* reads the XML file  */
//	int count = 0;
//	/* test if libxml works and if xml is parsed*/
//	LIBXML_TEST_VERSION;
//	xmlDocPtr doc;
//	doc = xmlReadFile(filename, NULL, 0);
//	if (doc == NULL)
//	{
//		fprintf(stderr, "Failed to parse %s\n", filename);
//		return false;
//	}
////finds the root of the xml tree
//	xmlNodePtr root_node = NULL, trace_node = NULL, line_node = NULL, bit_node = NULL;
//	root_node = xmlDocGetRootElement(doc);
//	//std::cout << "Root node: " << root_node->name <<std::endl;
//	
//	line_node =  root_node->children;
//	//reading lines
//	for ( ; line_node; line_node = line_node->next)
//	{		
//		if (line_node->type == XML_ELEMENT_NODE && !xmlStrcmp(line_node->name, BAD_CAST "TraceLine"))		
//		{		
//			float lid = xmlReadFloat(line_node,"ID",1);// default should never come for this case
//			int type = xmlReadInt(line_node,"Type",1);
//			int parent = xmlReadInt(line_node,"Parent",-1);
//
//			TraceLine * tline = new TraceLine();
//			tline->SetId(lid);
//			tline->SetType(type);
//			hash_load[lid]=reinterpret_cast<unsigned long long int>(tline);
//			if(parent!=-1)
//			{
//				TraceLine * tparent = reinterpret_cast<TraceLine*>(hash_load[parent]);
//				tline->SetParent(tparent);
//				tparent->GetBranchPointer()->push_back(tline);
//			}
//			else
//			{
//				trace_lines.push_back(tline);
//			}
//			//std::cout << "Line node: " << line_node->name << " " <<lid << std::endl;
//			int counter= 0, pid =0; double x, y, z;
//			bit_node = line_node->children;
//			
//			//reading points
//			for ( ;bit_node; bit_node = bit_node->next)
//			{
//				TraceBit tbit;
//				if ( !xmlStrcmp(bit_node->name, BAD_CAST "TraceBit") )	
//				{
//					pid = xmlReadInt(bit_node,"ID",0);
//					x = xmlReadFloat(bit_node,"X",0);
//					y = xmlReadFloat(bit_node,"Y",0);
//					z = xmlReadFloat(bit_node,"Z",0);
//					//printf("pid %d\n",pid);
//					tbit.x = x; tbit.y = y; tbit.z = z; tbit.id = pid;
//					tline->AddTraceBit(tbit);
//					counter++;			
//				}				
//			}
//			
//		}
//	}	
//	printf("Done\n");
//	return true;
//}


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

};
