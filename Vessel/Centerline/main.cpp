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

#include "vertex.h"
#include "edge.h"
#include "face.h"
#include "mesh.h"
#include <stdlib.h>
#include "radiosity.h"
#include <limits.h>
#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif 

//#pragma warning(disable : 4996)
//#pragma warning(disable : 4244)
//#pragma warning(disable : 4715)
//#pragma warning(disable : 4018)
//#pragma warning(disable : 4311)
//#pragma warning(disable : 4715)

#include "glCanvas.h"
#include "boundingbox.h"
#include <stdio.h>
#include <time.h>

#include "face.h"
#include <queue>

#ifndef M_PI
#define M_PI (2*acos(double(0)))
#endif

//#define M_PI (2*acos(0))
#define DEBUG printf
#define DEBUG_P 
#define epsilon 1e-3
#define epsilona 0.05
#define epsilonarea 0.5
//#include "marchtetra.h"
//#include "edge.h"
//#include "vertex.h"
//#include "face.h"
//#include "mesh.h"

using namespace std;

#include <vtksys/hash_map.hxx>

#define SKIP 2
#define MIN(a,b) (((a)>(b))?(b):(a))
#define MAX(a,b) (((a)>(b))?(a):(b))

//typedef hash_map<int ,Vertex *> type;
vtksys::hash_map< int ,Vertex * > vert;


vtksys::hash_map<long,double > hashcurv;
vtksys::hash_map<long,double > hashcurvnew;
vtksys::hash_map<long,bool > hash_decimate;

typedef vtksys::hash_multimap < int,Face *> map_type;
map_type hash_for_centerline;

Vec3f color1(0,0,0.5), color2(0,0,0.5), color3(0,0,0.5);

int *votes;
double lmin,lmedian,lmax;
int returned_count =0,non_returned_count=0;
unsigned char pmatrix[1030][1030][110];
//short int qmatrix[1030][1030][101];
bool debug_BFS= false;
double mincu, mediancu;
int numcu,totalcu;
Vec3f current_scaling;
double points[][3] = {
  {0, 0, 0},
  {1, 0, 0},
  {0, 1 ,0},
  {0, 0, 1},
  {0, 1, 1},
  {1, 0, 1},
  {1, 1, 0},
  {1, 1 ,1}
};

double total_min = 34324342,total_median = 0;
Face * total_min_face = NULL;

int tetra [][4]={
  {3,0,1,2},
  {3,5,7,1},
  {3,4,7,2},
  {3,7,2,1},
  {7,1,2,6},
  {0,5,1,6},
  {0,4,2,6},
  {0,4,6,5},
  {0,4,3,5},
  {7,4,6,5}
};     

struct twop{
	Vec3f start, destination;
};

vector<twop> debug_vector;
vector<twop> debug_vector1;
vector<twop> debug_vector2;
int find_max_around_current_point(Vec3f &);
void scale_vertices(double a, double b, double c);
void check_consistent(void);
void sleep(double delay1)
{
	clock_t start_time;
	start_time = clock();
	while((clock() - start_time) < delay1 * CLOCKS_PER_SEC)
	{
	}
}
double mratio = 0.5;
Mesh *m;        
int minx,miny,minz,maxx,maxy,maxz;

void CollectFaces(Vertex *have, Face *f, Array<Face*> *array) {
	if (array->Member(f)) return;
	if (have != (*f)[0] && have != (*f)[1] && have != (*f)[2]) return;
	array->Add(f);
	for (int i = 0; i < 4; i++) {
		Edge *ea = f->getEdge()->getOpposite();
		Edge *eb = f->getEdge()->getNext()->getOpposite();
		Edge *ec = f->getEdge()->getNext()->getNext()->getOpposite();
		// Edge *ed = f->getEdge()->getNext()->getNext()->getNext()->getOpposite();
		if (ea != NULL) CollectFaces(have,ea->getFace(),array);
		if (eb != NULL) CollectFaces(have,eb->getFace(),array);
		if (ec != NULL) CollectFaces(have,ec->getFace(),array);
		// if (ed != NULL) CollectFacesWithVertex(have,ed->getFace(),array);
	}
}

Vec3f ComputeN(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) {
	Vec3f v12 = p2;
	v12 -= p1;
	Vec3f v23 = p3;
	v23 -= p2;
	Vec3f normal;
	Vec3f::Cross3(normal,v12,v23);
	normal.Normalize();
	return normal;
}

float det2x2s(float a, float b,float c, float d) 
{
	return a * d - b * c;
}

float det3x3s(float a1,float a2,float a3,
			  float b1,float b2,float b3,
			  float c1,float c2,float c3) 
{
	return
		a1 * det2x2s( b2, b3, c2, c3 )
		- b1 * det2x2s( a2, a3, c2, c3 )
		+ c1 * det2x2s( a2, a3, b2, b3 );
}



Vec3f ComputeN(Face *f) {
	return ComputeN((*f)[0]->get(),(*f)[1]->get(),(*f)[2]->get());
}

Vec3f getCenter(Face * f)
{
	Edge *e = f->getEdge();
	Vec3f center = Vec3f(0,0,0);
	for(int co = 0; co <3; co++)
	{
		center = center + e->getVertex()->get();
		e = e->getNext();
	}
	center = center/3.0;
	return center;
}

float getLength(Edge *e){ 
	Vertex *v1,*v2;

	v1 = e->getVertex();
	v2 = e->getNext()->getVertex();
	Vec3f v= v1->get()-v2->get();
	return v.Length();
	return 0;
}


void BFS (Face *f, int depth, Array<Face*> *array)
{
	queue<Face*> q;
	queue<int> qindex;

	if(depth==0)
		return;

	Face *temp;
	Face * t;

	while(!q.empty())
		q.pop();
	while(!qindex.empty())
		qindex.pop();

	q.push(f);
	qindex.push(depth);
	while(!q.empty())
	{
		temp = q.front();
		int index = qindex.front();
		q.pop();
		qindex.pop();
		if(index==0)
		{
			continue;
		}
		Edge *e = temp->getEdge();

//		Edge * temp1 = NULL;
		for (int counter =0; counter < 3; counter ++)
		{


	//		Edge *e1 = e->getOpposite();
			if(e->getOpposite()!= NULL)
			{

				t = e->getOpposite()->getFace();
				if(!array->Member(t))
				{

					if(index>=1)
					{
						q.push(t);
						qindex.push(index-1);
					}
					array->Add(t);
				}
			}
			e=e->getNext();
		}
	}

}

void BFS_edge(Face *f, int depth, Array<Face*> *arrans,bool color)
{
	queue<Face*> q;
	queue<int> qindex;
	Array<Face*> array(70);
	if(depth==0)
		return;

	Face *temp;
	Face * t;

	while(!q.empty())
		q.pop();
	while(!qindex.empty())
		qindex.pop();

	q.push(f);
	qindex.push(depth);
	while(!q.empty())
	{
		temp = q.front();
		int index = qindex.front();
		q.pop();
		qindex.pop();
		if(index==0)
		{
			continue;
		}
		Edge *e = temp->getEdge();

//		Edge * temp1 = NULL;
		for (int counter =0; counter < 3; counter ++)
		{


	//		Edge *e1 = e->getOpposite();
			if(e->getOpposite()!= NULL)
			{

				t = e->getOpposite()->getFace();
				if(!array.Member(t))
				{

					if(index>=1)
					{
						q.push(t);
						qindex.push(index-1);
					}
					array.Add(t);
					if(index==1)
					{
						arrans->Add(t);
						if(color)
							t->setColor(Vec3f(0.0,0.0,1.0));
					}
				}
			}
			e=e->getNext();
		}
	}

}
void BFS1 (Face *f, Face *f1, int depth1, int depth2, Array<Face*> *array,bool color)
{
	queue<Face*> q;
	queue<int> qindex;

	if(depth2==0)
		return;

	Face *temp;
	Face * t;

	while(!q.empty())
		q.pop();
	while(!qindex.empty())
		qindex.pop();

	q.push(f1);
	qindex.push(depth2);
	while(!q.empty())
	{
		temp = q.front();
		int index = qindex.front();
		q.pop();
		qindex.pop();
		// printf("%d\n",index);
		if(index==0)
		{
			if(!q.empty())
				printf("wo no!\n");
			continue;
		}
		if(index ==3)
		{
			//	   printf("tempy");		
		}
		Edge *e = temp->getEdge();
		//printf("hi\n");
		//	   printf("%d\n",e);
		//Edge * temp1 = NULL;

		double dist = -10;
		Face * ftemp=NULL;
		for (int counter =0; counter < 3; counter ++)
		{

			//  printf("hi2\n");
			//Edge *e1 = e->getOpposite();
			// printf("hi-inter\n");
			if(e->getOpposite()!= NULL)
			{
				// printf("hi3\n");
				//fflush(stdin);
				t = e->getOpposite()->getFace();
				// printf("hi4\n");
				//if(!array->Member(t))
				if(!t->getMark())
				{
					//  if(t==NULL)
					//  	   printf(":OOOOOOOOOOOOO\n");
					Vec3f vtemp = getCenter(t)-getCenter(f);
					if(dist < vtemp.Length())
					{
						ftemp = t;
						dist = vtemp.Length();
					}

				}
			}
			//printf("test\n");
			e=e->getNext();
			//printf("test1\n");
		}
		if(ftemp!=NULL)
		{
			if(index>=1)
			{
				q.push(ftemp);
				qindex.push(index-1);
			}
			if(index <= depth2-depth1 )
			{
				array->Add(ftemp);
				if(color)
					ftemp->setColor(Vec3f(1.0,0,0));
				ftemp->setMark();
				//		ftemp->setColor(Vec3f(1.0,0.0,0.0));
			}
		}
	}

}

int BFS2(Face *f, Face *f1, int depth1, Array<Face*> *array)
{
	queue<Face*> q;
	queue<int> qindex;



	Face *temp;
	Face * t;

	while(!q.empty())
		q.pop();
	while(!qindex.empty())
		qindex.pop();

	q.push(f1);
	qindex.push(depth1);

	while(!q.empty())
	{
		temp = q.front();
		int index = qindex.front();
		q.pop();
		qindex.pop();
		// printf("%d\n",index);
		if(index==0)
		{
			break;
		}
		Edge *e = temp->getEdge();
		//printf("hi\n");
		//	   printf("%d\n",e);
//		Edge * temp1 = NULL;

		double dist = 1111110;
		Face * ftemp=NULL;
		for (int counter =0; counter < 3; counter ++)
		{

			//  printf("hi2\n");
	//		Edge *e1 = e->getOpposite();
			// printf("hi-inter\n");
			if(e->getOpposite()!= NULL)
			{
				// printf("hi3\n");
				//fflush(stdin);
				t = e->getOpposite()->getFace();
				// printf("hi4\n");
				//if(!array->Member(t))
				if(!t->getMark())
				{
					//  if(t==NULL)
					//  	   printf(":OOOOOOOOOOOOO\n");
					Vec3f vtemp = getCenter(t)-getCenter(f);
					if(dist > vtemp.Length())
					{
						ftemp = t;
						dist = vtemp.Length();
					}

				}
			}
			//printf("test\n");
			e=e->getNext();
			//printf("test1\n");
		}
		if(ftemp!=NULL)
		{
			if(index>=1)
			{
				q.push(ftemp);
				qindex.push(index-1);
			}
			array->Add(ftemp);
			ftemp->setMark();
			if(ftemp==f)
				return depth1-index+1;
			//		ftemp->setColor(Vec3f(1.0,0.0,0.0));
		}

	}
	return depth1;
}

bool BFS3(Face *f, Face *f1, int depth1, Array<Face*> *array)
{
	//	queue<Face*> q;
	//	queue<int> qindex;


	Array <Face*> arre(50);
	Face *temp;
	Face * t;

	/*	while(!q.empty())
	q.pop();
	while(!qindex.empty())
	qindex.pop();
	*/	
	//	q.push(f1);
	//	qindex.push(depth1);
	Vec3f current_direction = getCenter(f1)-getCenter(f);
	current_direction.Normalize();
	Vec3f new_direction = f->getNormal();
	new_direction.Normalize();
	Vec3f ellipse_normal ;
	Vec3f::Cross3(ellipse_normal,current_direction,new_direction);
	//	Vec3f sense = current_direction;

	Face* present_face=f1;
	int present_index = depth1;
	Vec3f ellipse_current;
	Vec3f vnormal;
	while(1)
	{
		temp = present_face;//q.front();
		//int index = qindex.front();
		//q.pop();
		//qindex.pop();
		// printf("%d\n",index);
		if(present_index==0)
		{
			break;
		}
		Edge *e = temp->getEdge();
		//printf("hi\n");
		//	   printf("%d\n",e);
		//Edge * temp1 = NULL;

		double dist = 2321312;
		Face * ftemp=NULL;
		arre.Clear();
		BFS(temp,4,&arre);
		vnormal = temp->getNormal();
		vnormal.Normalize();

		for (int counter =0; counter < arre.Count(); counter ++)
		{


			t = arre[counter];
			if(!t->getMark())
			{
				//  if(t==NULL)
				//  	   printf(":OOOOOOOOOOOOO\n");
				Vec3f vtemp = getCenter(t)-getCenter(temp);
				vtemp.Normalize();			
				if(dist > fabs(vtemp.Dot3(ellipse_normal))  ) 
				{
					Vec3f::Cross3(ellipse_current,vtemp,vnormal);
					if(ellipse_current.Dot3(ellipse_normal)>0)
					{
						ftemp = t;
						dist = fabs(vtemp.Dot3(ellipse_normal));
					}
				}
				//	printf("%d %d %lf \t",depth1,index,vtemp.Dot3(ellipse_normal));

			}

			//printf("test\n");
			e=e->getNext();
			//printf("test1\n");
		}
		if(ftemp!=NULL)
		{
			//		current_direction = getCenter(ftemp)-getCenter(temp);
			//		current_direction.Normalize();
			if(present_index>=1)
			{
				present_face = ftemp;
				present_index = present_index -1;
				//				q.push(ftemp);
				//			qindex.push(index-1);
			}
			array->Add(ftemp);
			ftemp->setMark();
			if(ftemp==f)
				return true;
			//		ftemp->setColor(Vec3f(1.0,0.0,0.0));
		}
		else
		{

			//	printf("I couldnt find ftemp\n");
			return false;
		}

	}
	return true;
}

bool BFS4(Face *f, Face *f1, int depth1, Array<Face*> *array,Array<Face*>*oldarray)
{
	//	queue<Face*> q;
	//	queue<int> qindex;


	Array <Face*> arre(50);
	Face *temp;
	Face * t;

	/*	while(!q.empty())
	q.pop();
	while(!qindex.empty())
	qindex.pop();
	*/	
	//	q.push(f1);
	//	qindex.push(depth1);
	Vec3f current_direction = getCenter(f1)-getCenter(f);
	current_direction.Normalize();
	Vec3f new_direction = f->getNormal();
	new_direction.Normalize();
	Vec3f ellipse_normal ;
	Vec3f::Cross3(ellipse_normal,current_direction,new_direction);
	//	Vec3f sense = current_direction;

	Face* present_face=f1;
	int present_index = depth1;
	Vec3f ellipse_current;
	Vec3f vnormal;
	while(1)
	{
		temp = present_face;//q.front();
		//int index = qindex.front();
		//q.pop();
		//qindex.pop();
		// printf("%d\n",index);
		if(present_index==0)
		{
			break;
		}
		Edge *e = temp->getEdge();
		//printf("hi\n");
		//	   printf("%d\n",e);
		//Edge * temp1 = NULL;

		double dist = 2321312;
		Face * ftemp=NULL;
		arre.Clear();
		BFS(temp,2,&arre);
		vnormal = temp->getNormal();
		vnormal.Normalize();

		for (int counter =0; counter < arre.Count(); counter ++)
		{


			t = arre[counter];

			if(!t->getMark())
			{
				//  if(t==NULL)
				//  	printf(":OOOOOOOOOOOOO\n");
				Vec3f vtemp = getCenter(t)-getCenter(temp);
				vtemp.Normalize();			
				if(dist > fabs(vtemp.Dot3(ellipse_normal))  ) 
				{
					Vec3f::Cross3(ellipse_current,vtemp,vnormal);
					if(ellipse_current.Dot3(ellipse_normal)>0)
					{
						if(oldarray->Member(ftemp))
						{
							ftemp = NULL;
							return false;
						}
						ftemp = t;
						dist = fabs(vtemp.Dot3(ellipse_normal));
					}
				}
				//	printf("%d %d %lf \t",depth1,index,vtemp.Dot3(ellipse_normal));

			}

			//printf("test\n");
			e=e->getNext();
			//printf("test1\n");
		}
		if(ftemp!=NULL)
		{


			//		current_direction = getCenter(ftemp)-getCenter(temp);
			//		current_direction.Normalize();
			if(present_index>=1)
			{
				present_face = ftemp;
				present_index = present_index -1;
				//				q.push(ftemp);
				//			qindex.push(index-1);
			}
			array->Add(ftemp);
			ftemp->setMark();
			if(ftemp==f)
				return true;
			//ftemp->setColor(Vec3f(1.0,0.0,0.0));
		}
		else
		{

			//	printf("I couldnt find ftemp\n");
			return false;
		}

	}
	return true;
}


int inline calc(Vec3f v);
void load_points( char * fname)
{

	minx = miny = minz = INT_MAX;
	maxx = maxy = maxz = INT_MIN;

	FILE *fp = fopen (fname,"r");
	int x,y,z;
	double z1,y1,x1;
	int count = 0;
	if(fp==NULL)
		printf("wo no! fp is null\n");

	while (!feof(fp))
	{
		if( fscanf(fp,"%lf %lf %lf %*f %*f %*f",&z1,&y1,&x1) == EOF )
      {
      cerr << "EOF encountered in fscanf" << endl;
      }
		//		if(!(x1 <300 &&  y1 <600 && y1>400 ))
		//				continue;

		z = (int)z1;y=(int)y1;x=(int)x1;
		count++;
		//	x = 1025-x;
		//	 y = 1025-y;
		//	z = 78-z;
		//	printf("%d %d %d\n",(x-1)/SKIP,(y-1)/SKIP,z-1);

		pmatrix[(x-1)/SKIP+1][(y-1)/SKIP+1][z-1+1]=1; // shrinking the output.
		minx = MIN(minx,(x-1)/SKIP);
		miny = MIN(miny,(y-1)/SKIP);
		minz = MIN(minz,z-1);
		maxx = MAX(maxx,(x-1)/SKIP+2);
		maxy = MAX(maxy,(y-1)/SKIP+2);
		maxz = MAX(maxz,z-1+2);
	}
	printf("Read %d points. maxx = %d, maxy =%d,maxz = %d\n",count,maxx, maxy,maxz);
	fclose(fp);
}

void load_points_without_normal( char * fname)
{

	minx = miny = minz = INT_MAX;
	maxx = maxy = maxz = INT_MIN;

	FILE *fp = fopen (fname,"r");
	int x,y,z;
	double z1,y1,x1;
	int l1;
	int count = 0;
	double l2;
	if(fp==NULL)
		printf("wo no! fp is null\n");

	while (!feof(fp))
	{

		if( fscanf(fp,"%lf %lf %lf %d %lf",&z1,&y1,&x1,&l1,&l2) == EOF )
      {
      cerr << "EOF encountered in fscanf" << endl;
      }
		//if(!(x1 <256 && y1<256))
		//	continue;
		//	if(!(x1 >400 && x1< 600  && y1>400 && y1 <600 ))
		//			continue;
		//if(!(x1 > 750 && x1 <950 && y1 > 75 && y1<450))
		//	continue;
		//if(!(x1 > 400 && x1 < 570 && y1 < 160))
		//	continue;
//		if(!(x1 > 470 && x1 < 600 && y1 < 370 && y1 > 230))
//			continue;
		// something
	//	if(!(x1 >830 && x1 < 970 && y1 >480 && y1 < 660 &&z1 >20 && z1 <30))
	//		continue;
		z = (int)z1;y=(int)y1;x=(int)x1;
		count++;
		//	x = 1025-x;
		//y = 1025-y;
		//	z = 78-z;
		//	printf("%d %d %d\n",(x-1)/SKIP,(y-1)/SKIP,z-1);

		pmatrix[(x-1)/SKIP+1][(y-1)/SKIP+1][z-1+1]=MAX(l1,pmatrix[(x-1)/SKIP+1][(y-1)/SKIP+1][z-1+1]);//((int(l2+0.5)==0)?l2-1.5:l2); // shrinking the output.
		minx = MIN(minx,(x-1)/SKIP);
		miny = MIN(miny,(y-1)/SKIP);
		minz = MIN(minz,z-1);
		maxx = MAX(maxx,(x-1)/SKIP+2);
		maxy = MAX(maxy,(y-1)/SKIP+2);
		maxz = MAX(maxz,z-1+2);
	}
	printf("Read %d points. maxx = %d, maxy =%d,maxz = %d\n",count,maxx, maxy,maxz);
	fclose(fp);
	//scanf("%*d");
}

//void read_from_file(char * filename)
//{
//	FILE * fp = fopen(filename,"r");
//	int ski = 2;
//	int linteger;
//	int maxlint = -1;
//	for (int countery =0; countery <1024; countery ++)
//		for (int counterx =0; counterx <1024; counterx ++)
//			for (int counterz =0; counterz < 58; counterz ++)
//			{
//				fscanf(fp,"%d",&linteger);
//				maxlint = MAX(maxlint,linteger);
//				qmatrix[counterx/ski+1][countery/ski+1][counterz+1]=linteger;
//			}
//			fclose(fp);
//			printf("maxlint %d\n",maxlint);
//}

//void load_points_from_qmatrix(int label)
//{
//	bool var;
//
//	for (int countery =0; countery <1024; countery ++)
//		for (int counterx =0; counterx <1024; counterx ++)
//			for (int counterz =0; counterz < 58; counterz ++)
//			{
//				var = (qmatrix[counterx][countery][counterz]==label);
//				pmatrix[counterx][countery][counterz] = var;
//				if(var)
//				{
//					minx = MIN(minx,(counterx));
//					miny = MIN(miny,(countery ));
//					minz = MIN(minz,counterz);
//					maxx = MAX(maxx,counterx+2);
//					maxy = MAX(maxy,countery+2);
//					maxz = MAX(maxz,counterz+2);
//				}
//			}
//}
void load_points_other( char * fname)
{
	minx = miny = minz = INT_MAX;
	maxx = maxy = maxz = INT_MIN;

	FILE *fp = fopen (fname,"r");


	//int count = 0;
	if(fp==NULL)
	{
		printf("wo no! fp is null\n");
		for(;;);
	}

	int linteger;


	//int ran =0;
	//char a=0;
	//bool flag = false;
	//	unsigned char c[2];



	for (int countery =0; countery <1024; countery ++)
		for (int counterx =0; counterx <1024; counterx ++)
			for (int counterz =0; counterz < 58; counterz ++)
			{
				//		int ret =fscanf(fp,"%c%c",&a,&b);
				//		fscanf(
				//		fread(c,1,1,fp);
				//		fread(c+1,1,1,fp);
				if( fscanf(fp,"%d",&linteger) == EOF )
          {
          cerr << "EOF encountered in fscanf" << endl;
          }
				//	if(ret!=2)
				//	{
				//		printf("yo! %d\n",rand());
				//	printf("%d %d %d\n",countery,counterx,counterz);
				//		for(;;);
				//	}
				//		if(linteger !=0)
				//			printf("hi\n");
				//		else
				//			printf("hi1\n");
				//		printf("%d",linteger);
				//
				//		linteger = c[0]+c[1]*256;
				//		if(linteger!=0 && linteger!=1)
				//			{
				//				printf("problem %d\n",linteger);
				//				for(;;);
				//			}
				//if(linteger)
				//		printf("%d ",linteger);
				//		fscanf(fp,"%c",&b);
				//	if(counterx < 400 && counterx > 100 && countery <700 && countery >400)
				int ski = 3;
				//if(linteger == 911 || linteger == 986 || linteger == 961 || linteger == 1010)
				if(linteger == 903 || linteger == 594)
				{
					pmatrix[(counterx/ski)+1][countery/ski+1][counterz+1]=(linteger!=0);
					minx = MIN(minx,(counterx)/ski);
					miny = MIN(miny,(countery/ski ));
					minz = MIN(minz,counterz);
					maxx = MAX(maxx,counterx/ski+2);
					maxy = MAX(maxy,countery/ski+2);
					maxz = MAX(maxz,counterz+2);
				}
				//	pmatrix[counterx/4][countery/4][counterz +1]=MAX(linteger!=0,pmatrix[counterx/4][countery/4][counterz+1]);
				//
				//		if(feof(fp))
				//		{
				////			printf("I already reached eof\n");
				//			flag = true;
				//		}
				//		else
				//		{
				//			if(flag)
				//				{
				//					printf("we're back!\n");
				//					for(;;);
				//				}
				//		}


			}

			//printf("I got %d points of label 150\n",ran);
			fclose(fp);

			//	while (!feof(fp))
			//	{
			//		fscanf(fp,"%lf %lf %lf",&z1,&y1,&x1);
			//		if(!(x1 >600 && x1 < 800 && y1 > 600 && y1 < 800 ))
			//				continue;
			//		z = (int)z1;y= (int)y1; x = (int)x1;
			//		count++;
			//		//	printf("%d %d %d\n",(x-1)/SKIP,(y-1)/SKIP,z-1);
			//
			//		pmatrix[(x-1)/SKIP+1][(y-1)/SKIP+1][z-1+1]=1; // shrinking the output.
			//		minx = MIN(minx,(x-1)/SKIP+1);
			//		miny = MIN(miny,(y-1)/SKIP+1);
			//		minz = MIN(minz,z-1+1);
			//		maxx = MAX(maxx,(x-1)/SKIP+2);
			//		maxy = MAX(maxy,(y-1)/SKIP+2);
			//		maxz = MAX(maxz,z-1+2);
			//	}
			printf("Read all the points points.\n");
			//scanf("%*d");
			//	minx = miny = minz = 1;
			//	maxx = maxy = 1025;
			//	maxz = 60;

			//	fclose(fp);
}



void marchtetra(void)
{
	vtksys::hash_map<int ,Vertex *>::iterator i;
	Vec3f a,b,c,d;
	Vec3f p,q,r;
	vector <Vec3f> arr;
	Vec3f res;
	Vec3f res1;
	Vec3f mid;
	Vertex *v1,*v2,*v3;
	vert.clear();
	for (int x = 0 ; x <maxx; x++  )
	{
		DEBUG("%d\n",x);
		//	 	scanf("%*d");
		for (int y = 0;  y<maxy ; y++ )
		{
			for(int z =0; z <maxz; z++)
			{
				//	 		sleep(0.05);
				//if(x>0)
				//DEBUG_P("%d %d %d\n",x,y,z);

				//for loop over 5 tetra hedras
				//add it to a list of faces with edges and vertices.
				for (int counter = 0; counter < 5; counter++)
				{
					if((x+y+z)%2)
					{
						a.Set(x+points[tetra[counter][0]][0],y+points[tetra[counter][0]][1],z+points[tetra[counter][0]][2]);
						b.Set(x+points[tetra[counter][1]][0],y+points[tetra[counter][1]][1],z+points[tetra[counter][1]][2]);
						c.Set(x+points[tetra[counter][2]][0],y+points[tetra[counter][2]][1],z+points[tetra[counter][2]][2]);
						d.Set(x+points[tetra[counter][3]][0],y+points[tetra[counter][3]][1],z+points[tetra[counter][3]][2]);
					}
					else
					{
						a.Set(x+points[tetra[counter+5][0]][0],y+points[tetra[counter+5][0]][1],z+points[tetra[counter+5][0]][2]);
						b.Set(x+points[tetra[counter+5][1]][0],y+points[tetra[counter+5][1]][1],z+points[tetra[counter+5][1]][2]);
						c.Set(x+points[tetra[counter+5][2]][0],y+points[tetra[counter+5][2]][1],z+points[tetra[counter+5][2]][2]);
						d.Set(x+points[tetra[counter+5][3]][0],y+points[tetra[counter+5][3]][1],z+points[tetra[counter+5][3]][2]);
					}
					int num = (pmatrix[(int)a.x()][(int)a.y()][(int)a.z()]!=0)+(pmatrix[(int)b.x()][(int)b.y()][(int)b.z()]!=0)+(pmatrix[(int)c.x()][(int)c.y()][(int)c.z()]!=0)+(pmatrix[(int)d.x()][(int)d.y()][(int)d.z()]!=0);

					//if(num!=0 && num !=4)
					//	printf("I did get some good nums\n");
					switch(num)
					{
					case 0:
					case 4:
						break;
					case 3:
						//DEBUG_P("case 3\n");
						if(pmatrix[(int)a.x()][(int)a.y()][(int)a.z()]==0)
						{
							p=a*(1-mratio)+b*mratio;
							q=a*(1-mratio)+c*mratio;
							r=a*(1-mratio)+d*mratio;
							mid = a;

						}
						else if (pmatrix[(int)b.x()][(int)b.y()][(int)b.z()]==0)
						{
							p=a*mratio+b*(1-mratio);
							q=b*(1-mratio)+c*mratio;
							r=b*(1-mratio)+d*mratio;
							mid = b;
						}
						else if (pmatrix[(int)c.x()][(int)c.y()][(int)c.z()]==0)
						{
							p=a*mratio+c*(1-mratio);
							q=b*mratio+c*(1-mratio);
							r=c*(1-mratio)+d*mratio;
							mid = c;
						}
						else
						{
							p=a*mratio+d*(1-mratio);
							q=b*mratio+d*(1-mratio);
							r=c*mratio+d*(1-mratio);
							mid = d;
						}

						Vec3f::Cross3(res,q-p,r-q);

						Vec3f::Mult(res1,res,(p+q+r)-3*mid);
						if(res1.x()+res1.y()+res1.z()>0)
						{
							Vec3f swap = p;
							p = q;
							q = swap;
						}
						i = vert.find(calc(p));
						//Vertex *v1;
						if(i !=vert.end())
						{
							v1 = i->second;  
						}
						else
						{
							v1 = m->addVertex(p);
							vert[calc(p)]=v1;
						}
						i = vert.find(calc(q));
						// Vertex *v2;
						if(i !=vert.end())
						{
							v2 = i->second;  
						}
						else
						{
							v2 = m->addVertex(q);
							vert[calc(q)]=v2;
						}
						i = vert.find(calc(r));
						// Vertex *v3;
						if(i !=vert.end())
						{
							v3 = i->second;  
						}
						else
						{
							v3 = m->addVertex(r);
							vert[calc(r)]=v3;
						}
						m->addFace(v1,v2,v3,Vec3f(1,1,1),color1);
						break;

					case 1:
						//DEBUG_P("case 1\n");
						if(pmatrix[(int)a.x()][(int)a.y()][(int)a.z()]!=0)
						{
							//DEBUG_P("TT-a");
							p=a*(mratio)+b*(1-mratio);
							q=a*(mratio)+c*(1-mratio);
							r=a*(mratio)+d*(1-mratio);
							mid = a;

						}
						else if (pmatrix[(int)b.x()][(int)b.y()][(int)b.z()]!=0)
						{
							//DEBUG_P("TT-b");
							p=a*(1-mratio)+b*(mratio);
							q=b*(mratio)+c*(1-mratio);
							r=b*(mratio)+d*(1-mratio);
							mid = b;
						}
						else if (pmatrix[(int)c.x()][(int)c.y()][(int)c.z()]!=0)
						{
							//DEBUG_P("TT-c");
							p=a*(1-mratio)+c*(mratio);
							q=b*(1-mratio)+c*(mratio);
							r=c*(mratio)+d*(1-mratio);
							mid = c;
						}
						else
						{
							//DEBUG_P("TT-d");
							p=a*(1-mratio)+d*(mratio);
							q=b*(1-mratio)+d*(mratio);
							r=c*(1-mratio)+d*(mratio);
							mid = d;
						}

						Vec3f::Cross3(res,q-p,r-q);

						Vec3f::Mult(res1,res,(p+q+r)-3*mid);
						if(res1.x()+res1.y()+res1.z()<0)
						{
							//DEBUG_P("I did come here\n");
							Vec3f swap = p;
							p = q;
							q = swap;
						}
						i = vert.find(calc(p));

						if(i !=vert.end())
						{
							v1 = i->second;  
						}
						else
						{
							v1 = m->addVertex(p);
							vert[calc(p)]=v1;
						}
						i = vert.find(calc(q));
						//Vertex *v2;
						if(i !=vert.end())
						{
							v2 = i->second;  
						}
						else
						{
							v2 = m->addVertex(q);
							vert[calc(q)]=v2;
						}
						i = vert.find(calc(r));
						//Vertex *v3;
						if(i !=vert.end())
						{
							//DEBUG_P("its because I came here \n");
							v3 = i->second;  
						}
						else
						{
							v3 = m->addVertex(r);
							vert[calc(r)]=v3;
						}
						/*
						if(calc(v1->get())==calc(v2->get())||calc(v1->get())==calc(v3->get())||calc(v2->get())==calc(v3->get()))
						{
						printf("two are same %d %d %d\n",calc(v1->get()),calc(v2->get()),calc(v3->get()));
						(printf("positions"),v1->print(),v2->print(),v3->print());
						(printf("a b c d"),a.Write(),b.Write(),c.Write(),d.Write());
						printf("check");
						(c*(1-mratio)+d*(mratio)).Write();
						printf("check");
						(p.Write(),q.Write(),r.Write());
						printf("maxx %d maxy %d maxz %d\n",maxx,maxy,maxz);
						}
						*/
						m->addFace(v1,v2,v3,Vec3f(1,1,1),color2);
						break;
					case 2:
						//DEBUG_P("case 2");
						arr.clear();

						if(pmatrix[(int)a.x()][(int)a.y()][(int)a.z()]!=0 && pmatrix[(int)b.x()][(int)b.y()][(int)b.z()]!=0)
						{
							arr.push_back(a*mratio+c*(1-mratio));
							arr.push_back(a*mratio+d*(1-mratio));
							arr.push_back(b*mratio+d*(1-mratio));
							arr.push_back(b*mratio+c*(1-mratio));
							mid = (a+b)/2.0;
						}
						if(pmatrix[(int)a.x()][(int)a.y()][(int)a.z()]!=0 && pmatrix[(int)c.x()][(int)c.y()][(int)c.z()]!=0)
						{
							arr.push_back(a*mratio+b*(1-mratio));
							arr.push_back(a*mratio+d*(1-mratio));
							arr.push_back(c*mratio+d*(1-mratio));
							arr.push_back(c*mratio+b*(1-mratio));
							mid = (a+c)/2.0;
						}
						if(pmatrix[(int)a.x()][(int)a.y()][(int)a.z()]!=0 && pmatrix[(int)d.x()][(int)d.y()][(int)d.z()]!=0)
						{
							arr.push_back(a*mratio+c*(1-mratio));
							arr.push_back(a*mratio+b*(1-mratio));
							arr.push_back(d*mratio+b*(1-mratio));
							arr.push_back(d*mratio+c*(1-mratio));
							mid = (a+d)/2.0;
						}
						if(pmatrix[(int)c.x()][(int)c.y()][(int)c.z()]!=0 && pmatrix[(int)b.x()][(int)b.y()][(int)b.z()]!=0)
						{
							arr.push_back(c*mratio+a*(1-mratio));
							arr.push_back(c*mratio+d*(1-mratio));
							arr.push_back(b*mratio+d*(1-mratio));
							arr.push_back(b*mratio+a*(1-mratio));
							mid = (b+c)/2.0;
						}
						if(pmatrix[(int)d.x()][(int)d.y()][(int)d.z()]!=0 && pmatrix[(int)b.x()][(int)b.y()][(int)b.z()]!=0)
						{
							arr.push_back(b*mratio+a*(1-mratio));
							arr.push_back(b*mratio+c*(1-mratio));
							arr.push_back(d*mratio+c*(1-mratio));
							arr.push_back(d*mratio+a*(1-mratio));
							mid = (b+d)/2.0;
						}
						if(pmatrix[(int)c.x()][(int)c.y()][(int)c.z()]!=0 && pmatrix[(int)d.x()][(int)d.y()][(int)d.z()]!=0)
						{
							arr.push_back(d*mratio+a*(1-mratio));
							arr.push_back(d*mratio+b*(1-mratio));
							arr.push_back(c*mratio+b*(1-mratio));
							arr.push_back(c*mratio+a*(1-mratio));
							mid = (c+d)/2.0;
						}
						vector <Vertex*> ps;
						Vec3f::Cross3(res,arr[1]-arr[0],arr[2]-arr[1]);
						Vec3f::Mult(res1,res,(arr[0]+arr[1]+arr[2]+arr[3])-4*mid);
						if(res1.x()+res1.y()+res1.z()<0)
						{
							Vec3f swap = arr[1];
							arr[1]=arr[3];
							arr[3]=swap;
						}
						for (int co = 0; co<4; co++)
						{
							i = vert.find(calc(arr[co]));
							Vertex *v;
							if(i !=vert.end())
							{
								v = i->second;  
							}
							else
							{
								v = m->addVertex(arr[co]);
								vert[calc(arr[co])]=v;
							}
							ps.push_back(v);
						}

						//if(unique(ps.begin(),ps.end())!=ps.end())
						//{printf("bg!\n");for(;;);}
						//if(unique(arr.begin(),arr.end())!=arr.end())
						//{printf("bg!\n");for(;;);}
						m->addFace(ps[0],ps[1],ps[2],Vec3f(1,1,1),color3);
						m->addFace(ps[0],ps[2],ps[3],Vec3f(1,1,1),color3);
						break;
					}
				}
			}
		}
	}
	printf("Finished - Marchtetra\n");
	// start from a single face... fix its orientation -> do a bfs aross the faces and fix the orientations -> grow the mesh
}

void Shift(Face *f,Vec3f center)
{
	Edge *e = f->getEdge();
	Vertex *v;
	Vec3f oldcenter = getCenter(f);
	for (int co = 0; co <3; co++)
	{
		v = e->getVertex();
		v->set(v->get()+center-oldcenter);
	}
	Vec3f temp = getCenter(f)-center;
	if(temp.Length()>0.01)
	{
		printf("praaaaaalam\n");
		for(;;);

	}
}
int inline calc(Vec3f v)
{	
	Vec3f vc = v;
	double x1=vc.x(),y1=vc.y(),z1=vc.z();
	if(fabs(vc.x()-int(vc.x()+0.5))>1e-4)
	{
		x1 = int(vc.x())+0.5;
	}
	if(fabs(vc.y()-int(vc.y()+0.5))>1e-4)
	{
		y1 = int(vc.y())+0.5;
	}
	if(fabs(vc.z()-int(vc.z()+0.5))>1e-4)
	{
		z1 = int(vc.z())+0.5;
	}
	//printf(" I got %lf %lf %lf x1 %lf y1 %lf z1 %lf\n",v.x(),v.y(),v.z(),x1,y1,z1);
	//printf(" I returned %d\n",int(int(2*x1+0.5)+int(2*y1+0.5)*(2*maxx+3)+int(2*z1+0.5)*(2*maxx+3)*(2*maxy+6)));
	return int(int(2*x1+0.5)+int(2*y1+0.5)*(2*maxx+3)+int(2*z1+0.5)*(2*maxx+3)*(2*maxy+6));
}

int inline calcint(Vec3f &v)
{
	Vec3f new1 = v;
	new1.Set(int(new1.x()+0.5),int(new1.y()+0.5),int(new1.z()+0.5));
	return int(new1.x()+new1.y()*(1024+1)+new1.z()*(1024+1)*(1024+1));
}

double FindCurvature(Face *f1,Face *f2)
{
	static int infcount = 0;
	Vec3f c1,c2;
	c1 = getCenter(f1);
	c2 = getCenter(f2);
	Vec3f n1,n2;
	n1 = ComputeN(f1);
	n2 = ComputeN(f2);
	Vec3f result;
	Vec3f::Cross3(result,n1,n2);
	Vec3f diff = c1-c2;
	double sign;

	n1 = 0.5*diff.Length()*n1;
	n2 = 0.5*diff.Length()*n2;

	Vec3f temp = c1 + n1 - ( c2 + n2);
	if(diff.Length() < temp.Length())
	{		   sign = 1;}
	else
	{
		//printf("I came here\n");
		sign = -1;
	}

	if(diff.Length() > epsilon)
		return sign*result.Length()/diff.Length();
	infcount ++;
	printf("yes (%lf %lf %lf)\t(%lf %lf %lf) %d!\n",c1.x(),c1.y(),c1.z(),c2.x(),c2.y(),c2.z(),infcount);
	return sign*1000;
}

Vec3f getColorCode(double c,double median,double minc, double maxc)
{
	Vec3f ans;


	if(c < median)
	{
		c= 0.5*( c- minc)/(median-minc);
		c = MAX(c,0);
		ans.Set(1-c*2,c*2,0);
	}
	else
	{
		c = 0.5 + 0.5*(c-median )/(maxc-median);
		c = MIN(c,1);
		ans.Set(0,2-2*c,2*c-1);
	}

	return ans;
}

double getDet(Face * f, Vertex * v)
{

	Edge * test = f->getEdge();
	while(test->getVertex()!=v)
	{
		test=test->getNext();
	}
	Vec3f v0 = test->getVertex()->get();
	Vec3f v1 = test->getNext()->getVertex()->get();
	Vec3f v2 = test->getNext()->getNext()->getVertex()->get();
	Vec3f res;
	Vec3f::Cross3(res,v1,v2);
	return res.Dot3(v0);
	//	return v0.x()*(v1.y()*v2.z()-v1.z()*v2.y())-v1.x()*(v0.y()*v2.z()-v0.z()*v2.y())+v2.x()*(v0.y()*v1.z()-v1.y()*v0.z());
}
void SmoothSurface(bool selective=false)
{

	Bag<Face *>*bf = m->getFaces();
	Iterator<Face*>* iter = bf->StartIteration();
	Array<Face*> arr(50);

	int rand_count=0;
	int solved_count = 0;
	int didnt_solve_count =0;
	if(iter == NULL)
	{
		printf("what the ... ?\n");
		for(;;);
	}
	while(1)
	{

		arr.Clear();
		Face *fstart = iter->GetNext();//bf->ChooseRandom();
		//	printf("I've come here\n");
		if(fstart==NULL)
			break;
		if(selective && !fstart->getMark())
			continue;
		//	fstart->setColor(Vec3f(1.0,1.0,1.0));
		//BFS(fstart,1,&arr);
		Edge *e = fstart->getEdge();

		Vertex *v1[3];
		v1[0] = e->getVertex();
		v1[1] = e->getNext()->getVertex();
		v1[2] = e->getNext()->getNext()->getVertex();

		for (int co = 0;co<3; co++)
		{
			arr.Clear();
			//	scanf("%*d");
			CollectFaces(v1[co],fstart, &arr);
			//BFS(fstart,4,&arr);
			//DEBUG("ha!\n");
			Vec3f sum = Vec3f(0,0,0);
			int sum_counter = 0;
			Vec3f Vnormal = Vec3f(0.0,0.0,0.0);
			//double Det_sum = 0;
			for( int counter =0; counter < arr.Count(); counter ++)
			{
				Vec3f center= getCenter(arr[counter]);
				sum = sum + center;
				sum_counter++;
				//				Vnormal += (arr[counter]->getNormal());
				//				Det_sum += getDet(arr[counter],v1[co]);
			}
			rand_count++;
			//	Vec3f Vm = sum/sum_counter;
			//Vnormal /=sum_counter;
			//		if(Vnormal.Length()<1e-1)
			//			{
			//				v1[co]->set(Vm);
			//				continue;
			//			}
			//DEBUGGING MODE

			v1[co]->set(sum/sum_counter);
			continue;
			// END OF DEBUGGING MODE
			//	double K = (Vnormal.Dot3(Vm)-Det_sum)/(Vnormal.Dot3(Vnormal));
			//	double K2 = (Vnormal.Dot3(Vm)+Det_sum)/(Vnormal.Dot3(Vnormal));

			//	if(fabs(K)>fabs(K2))
			//		K = K2;

			//		v1[co]->set(Vm-K*Vnormal);
			//		continue;
			//printf("%lf\n",K);

			//Shift(fstart,sum/sum_counter);
			//v1[co]->set(sum/sum_counter);
		}

		//if(rand_count %10 ==0)
		//printf("Done %d\n",rand_count);
	}
	bf->EndIteration(iter);
	printf("Solved count = %d Didnt Solve count = %d\n",solved_count, didnt_solve_count);
}

void SmoothSurface_volume_preserve()
{
	double total_abs_error=0;

	Bag<Face *>*bf = m->getFaces();
	Iterator<Face*>* iter = bf->StartIteration();
	Array<Face*> arr(50);

	int rand_count=0;
	//int solved_count = 0;
	//int didnt_solve_count =0;
	while(1)
	{

		arr.Clear();
		Face *fstart = iter->GetNext();//bf->ChooseRandom();

		if(fstart==NULL)
			break;

		fstart->setColor(Vec3f(0.0,1.0,0.0));
		//BFS(fstart,1,&arr);
		Edge *e = fstart->getEdge();

		Vertex *v1[3];
		v1[0] = e->getVertex();
		v1[1] = e->getNext()->getVertex();
		v1[2] = e->getNext()->getNext()->getVertex();

		for (int co = 0;co<3; co++)
		{
			arr.Clear();
			//	scanf("%*d");
			CollectFaces(v1[co],fstart, &arr);
			//BFS(fstart,4,&arr);
			//DEBUG("ha!\n");
			Vec3f sum = Vec3f(0,0,0);
			int sum_counter = 0;
			Vec3f Vnormal = Vec3f(0.0,0.0,0.0);
			double Det_sum = 0;
			vector <Vec3f> vec;
			for( int counter =0; counter < arr.Count(); counter ++)
			{
				Vec3f center= getCenter(arr[counter]);
				sum = sum + center;
				sum_counter++;
				Vnormal = Vnormal -(arr[counter]->getNormal());
				vec.push_back(Vnormal);
				Det_sum += getDet(arr[counter],v1[co]);
			}

			rand_count++;
			Vec3f Vm = sum/sum_counter;
			//Vnormal /=sum_counter;
			//		if(Vnormal.Length()<1e-1)
			//			{
			//				v1[co]->set(Vm);
			//				continue;
			//			}
			//DEBUGGING MODE
			//	v1[co]->set(Vm);
			//	continue;
			// END OF DEBUGGING MODE
			double K = (Vnormal.Dot3(Vm)-Det_sum)/(Vnormal.Dot3(Vnormal));
			double K2 = (+Vnormal.Dot3(Vm)-Det_sum)/(Vnormal.Dot3(Vnormal));

			//		if(fabs(-(v1[co]->get().Dot3(Vnormal))+Det_sum)<epsilon)
			//			{
			//				printf("K supposedly works %lf %lf\n",K,K2);
			//			}
			//			else if(fabs(v1[co]->get().Dot3(Vnormal)+Det_sum) < 1e-1)
			//			{
			//				printf("K2 supposedly works %lf %lf\n",K,K2);
			//			}
			//			else
			//			{
			//				printf("%lf\t",fabs(-(v1[co]->get().Dot3(Vnormal))+Det_sum));
			//				printf("%lf\t",fabs(v1[co]->get().Dot3(Vnormal)+Det_sum));
			//				for(int t = 0; t < vec.size(); t++)
			//				{
			//					for(int t1 = 0; t1 < vec.size(); t1++)
			//					{
			//						if(vec[t].Length()>epsilon && vec[t1].Length()>epsilon)
			//							printf("Angle %lf\t",acos(vec[t].Dot3(vec[t1])/vec[t].Length()/vec[t1].Length()));
			//					}
			//				}
			//			 	printf("  Nothing Works..... worsht! :(\n");
			//			 	scanf("%*d");
			//			}

			//	if(fabs(K)>fabs(K2))
			K = K2;
			Vec3f diff = v1[co]->get()-(Vm-K*Vnormal);
			total_abs_error += diff.Length();
			v1[co]->set(Vm-K*Vnormal);
			continue;
			//printf("%lf\n",K);

			//Shift(fstart,sum/sum_counter);
			//v1[co]->set(sum/sum_counter);
		}

		//if(rand_count %10 ==0)
		//printf("Done %d\n",rand_count);
	}
	bf->EndIteration(iter);
	//	printf("Solved count = %d Didnt Solve count = %d\n",solved_count, didnt_solve_count);
	printf("total abs error %lf\n",total_abs_error);
}

void SmoothSurface_complex( float w) // two point smoothing
{
	Bag <Edge *> * be = m->getEdges();
	Iterator<Edge *> * iter = be->StartIteration();
	while(Edge * e = iter->GetNext())
	{
		if(e->getOpposite()==NULL)
		{
			printf("This should not happen .. :-? but I took care of it! :) \n");
			e->getFace()->setColor(Vec3f(1,0,0));
			continue;
		}
		//smooth_edge(e);
		Edge * eiter = e;
		vector <Vec3f> v1;
		vector <Vec3f> v2;
		v1.clear();//defensive programming :P
		v2.clear();
		Vec3f x1 = (*e)[0]->get();
		Vec3f x2 = (*e)[1]->get();
		do
		{
			v1.push_back(((*eiter)[1])->get());
			if(eiter->getPrev()->getOpposite()==NULL)
			{
				// we probably shouldnt meddle with this boundary thingies
				// continue with anothe edge;
				break;
			}
			eiter=eiter->getPrev()->getOpposite();
		}while(eiter!=e);
		if(eiter!=e)
			continue;
		//set e to e->getOpposite() and do the same thing again
		e = e->getOpposite();// cant be null
		eiter = e;
		do
		{
			v2.push_back(((*eiter)[1])->get());
			if(eiter->getPrev()->getOpposite()==NULL)
			{
				// we probably shouldnt meddle with this boundary thingies
				// continue with anothe edge;
				break;
			}
			eiter=eiter->getPrev()->getOpposite();
		}while(eiter!=e);			
		if(eiter!=e)
			continue;

		//if(v1.size() < 5 && v2.size()<5)
		//continue;
		Vec3f A1=Vec3f(0,0,0),A2=Vec3f(0,0,0);
		Vec3f res;

		//find A1,A2
		for(unsigned int counter = 0; counter < v1.size(); counter ++)
		{
			Vec3f::Cross3(res,v1[counter]-x1,v1[(counter+1)%v1.size()]-x1);
			A1 = A1 + res;
		}
		for(unsigned int counter = 0; counter < v2.size(); counter ++)
		{
			Vec3f::Cross3(res,v2[counter]-x2,v2[(counter+1)%v2.size()]-x2);
			A2 = A2 + res;
		}

		//find v
		//		Vec3f v  = v1[1]-v1[v1.size()-1];
		Vec3f v  = v2[v2.size()-1]-v2[1];
		if((v1[0]-x2).Length()>1e-3 || (v2[0]-x1).Length()>1e-3)
		{
			printf("wraaang\n");
		}
		//find smoothed vertex
		Vec3f x1s = Vec3f(0,0,0);
		Vec3f x2s = Vec3f(0,0,0);

		for(unsigned int counter =1; counter < v1.size(); counter++)
			x1s += v1[counter];
		x1s = x1s * v2.size();
		for(unsigned int counter = 1; counter < v2.size(); counter++)
			x1s += v2[counter];

		x1s /= float(v1.size()*v2.size()-1);

		for(unsigned int counter = 1; counter < v2.size(); counter++)
			x2s += v2[counter];
		x2s = (x1s + x2s)/v2.size();

		Vec3f dx1s = w*(x1s - x1);
		Vec3f dx2s = w*(x2s - x2);

		Vec3f temp;
		Vec3f::Cross3(temp,v,dx1s-dx2s);
		Vec3f A = A1 + A2 + temp;
		double length = ((*e)[1]->get()-(*e)[0]->get()).Length();
		if(A.Length() > 0.1)
		{
			//	printf("%lf\t",A.Length());
			Vec3f normal = A/A.Length();
			Vec3f::Cross3(temp,v,dx1s);
			double h = -(dx1s.Dot3(A1) + dx2s.Dot3(A2) + dx2s.Dot3(temp))/A.Length();
			if((dx1s+h*normal).Length()<length&&(dx2s+h*normal).Length()<length)
			{
				x1 = x1 + dx1s + h*normal;
				x2 = x2 + dx2s + h*normal;
				e = e->getOpposite();
				(*e)[0]->set(x1);
				(*e)[1]->set(x2);
			}
		}
	}
	be->EndIteration(iter);
}

void FillCurvatures(int num)
{
	Bag<Face *>*bf = m->getFaces();
	Iterator<Face*>* iter ;
	Array<Face*> arr(50);

	iter = bf->StartIteration();
	double maxc =-1,minc=242123424;
	vector <double> curvs;
	curvs.clear();
	int trys = 0;
	while(Face *f = iter->GetNext())
	{
		arr.Clear();
		BFS(f,num,&arr);
		//	printf("size %d\n",arr.Count());
		if(trys%1000==0)
			printf("%d\n",trys);
		trys ++;
		//	if(arr.Count()!=3)
		//  	printf("something wrong!!!");

		double sum_curvature=0;
		int sum_count=0;
		for (int co = 0; co <arr.Count();co++)
		{
			if(f!=arr[co])
			{
				sum_curvature += FindCurvature(f,arr[co]);
				sum_count++;
			}
		}
		if(sum_count==0)
			continue;
		sum_curvature/=sum_count;
		curvs.push_back(sum_curvature);
		maxc=MAX(maxc,sum_curvature);
		minc=MIN(minc,sum_curvature);
		hashcurv[reinterpret_cast<long>(f)]=sum_curvature;
		//printf("%lf\t",sum_curvature/sum_count);
	}
	sort(curvs.begin(),curvs.end());
	double median = curvs[int(curvs.size()/2)];
	maxc = curvs[int(0.85*curvs.size())];minc = curvs[int(0.1*curvs.size())];
	printf("%lf %lf %lf\n",minc,median,maxc);
	bf->EndIteration(iter);

	iter = bf->StartIteration();
	while(Face *f = iter->GetNext())
	{
		Vec3f color = getColorCode(hashcurv[(long)f],median,minc,maxc);
		f->setColor(color);
		// 		   f->setColor(get
		//printf("%lf\t",sum_curvature/sum_count);
	}

	bf->EndIteration(iter);
	DEBUG("ENDED");
}

void FillCurvatures_rough(int num1,int num2,int num3)
{
	hashcurv.clear();
	hashcurvnew.clear();
	mincu = 123434;
	mediancu = 1023214;
	numcu = 0;
	totalcu = 0;

	Bag<Face *>*bf = m->getFaces();

	Iterator<Face*>* iter =bf->StartIteration();
	while(Face* f2 = iter->GetNext())
	{
		f2->unMark();
	}
	bf->EndIteration(iter);

	Array<Face*> arr(50);
	Array<Face*> arr2(50);
	iter = bf->StartIteration();

	double maxc =-1,minc=242123424;
	vector <double> curvs;
	curvs.clear();
	int seedpoints  = 10;
	int trys = 0;

	while(1)
	{
		Face *f = iter->GetNext();
		if(f==NULL)
			break;
		if(trys%1000==0)
			printf("%d\n",trys);
		//			if(trys > 20)
		//			break;
		trys ++;
		arr.Clear();

		if(trys>40 && trys <60)
			BFS_edge(f,num3,&arr,true&&debug_BFS);
		else
			BFS_edge(f,num3,&arr,false&&debug_BFS);	 		   

		double curvature =0;
		bool flag = false;
		double sum_curvature;
		int sum_count;

		for(int counter =0; counter < arr.Count(); counter++)
		{
			arr2.Clear();
			sum_curvature=0;
			sum_count = 0;
			//	if(trys>40 && trys <60)
			//		BFS1(f,arr[counter],num1,num2,&arr2,true);
			//	else
			BFS1(f,arr[counter],num1,num2,&arr2,false);
			for(int counter2 = 0; counter2 < arr2.Count(); counter2++)
			{
				arr2[counter2]->unMark();
				Vec3f getdiff = getCenter(f)-getCenter(arr2[counter2]);
				if(getdiff.Length()>epsilon)
				{
					sum_curvature += FindCurvature(f,arr2[counter2]);
					sum_count ++;
				}
			}
			if(sum_count==0)
				continue;

			//CHANGED
			if(flag)
			{
				curvature = ( fabs(curvature)>fabs(sum_curvature/sum_count))?curvature:sum_curvature/sum_count;
			}
			else
			{
				curvature = sum_curvature/sum_count;
				flag = true;
			}
			//					   curvature = MAX(curvature,sum_curvature/sum_count);
		}
		//   if(curvature < 0)
		//   				printf("I got -ve curvatures too!\n");
		curvs.push_back(curvature);
		maxc=MAX(maxc,curvature);
		minc=MIN(minc,curvature);
		hashcurv[(long)f]=curvature;
		seedpoints--;

		//   if(seedpoints <0)
		//  break;
	}
	sort(curvs.begin(),curvs.end());
	double median = curvs[int(curvs.size()/2)];
	mediancu = median;

	total_median = median;
	maxc = curvs[int(0.85*curvs.size())];minc = curvs[int(0)];
	minc = -0.7; maxc = 0.7; 
	median = 0.2;
	printf("%lf %lf %lf\n",minc,median,maxc);
	bf->EndIteration(iter);
	if(debug_BFS)
		return;
	iter = bf->StartIteration();
	int try1 = 0;
	while(Face *f = iter->GetNext())
	{
		try1++;
		if(try1%1000==0)
			printf("%d\r",try1);
		arr.Clear();
		BFS(f,6,&arr);
		vector <double> values;

		for(int counter =0; counter < arr.Count(); counter++)
		{
			values.push_back(hashcurv[(long)arr[counter]]);
		}
		sort(values.begin(),values.end());
		double median1 = values[values.size()/2];
		hashcurvnew[(long)f]=median1;
		if(mincu > median1)
			mincu = median1;
		if(median1 <=0)
			numcu++;

		totalcu++;;
		Vec3f color = getColorCode(median1,median,minc,maxc);
		if(median1<median)
		{
			//f->setColor(Vec3f(0.0,0.0,0.0));
			//	f->setEmit(Vec3f(0.0,0.0,0.0));
		}
		if(total_min > median1)
		{
			total_min = median1;
			total_min_face = f;
		}
		f->setColor(color);
		//		f->setColor(Vec3f(0,1,0));
		//	continue;
		//		if(!debug_BFS)
		//			if(median1<=-10000000)
		//				f->setColor(Vec3f(1.0,0.0,0.0));
		//			else
		//			{
		//				if(median1<-0.1)
		//					f->setColor(Vec3f(1.0,0.0,0.0));
		//				else if(median1 >0.1)
		//					f->setColor(Vec3f(0,0,1.0f));
		//				else
		//					f->setColor(Vec3f(0,1,0));
		////				f->setColor(color);
		//			}
		// 		   f->setColor(get
		//printf("%lf\t",sum_curvature/sum_count);

	}

	bf->EndIteration(iter);

}

void rank70Filter()
{
	printf("Started 70%% rank filter\n");
	int mat[27][3];
	int matcount=0;
	vector <int> a;
	for(int counterx =-1; counterx < 2; counterx ++)
		for(int countery = -1; countery < 2; countery ++ )
			for(int counterz = -1; counterz <2 ; counterz ++)
			{
				mat[matcount][0]=counterx;
				mat[matcount][1]=countery;
				mat[matcount][2]=counterz;
				matcount ++;
			}

			for(int counterx = 1;counterx <maxx-1; counterx ++)
				for(int countery = 1; countery < maxy-1; countery ++)
					for(int counterz = 0; counterz < maxz; counterz++)
					{

						int sum = 0;
						for(int counterw = 0; counterw < 27; counterw++)
						{
							if(counterw!=13)
							{
								int temp1 = MAX(0,mat[counterw][2]+counterz);
								temp1 = (int)MIN(temp1,maxz-1);
								//a.push_back(pmatrix[counterx][countery][temp1]);
								sum = sum + pmatrix[counterx+mat[counterw][0]][countery+mat[counterw][1]][temp1];
							}
						}
						if(sum > 7)
							pmatrix[counterx][countery][counterz]=1;
						else
							pmatrix[counterx][countery][counterz]=0;
					}

					printf("Finished 70%% rank filter\n");
}

void removeInnerSurface()
{

	if(total_min >0)
		return;
	total_median = MAX(0,total_median);

	// Start with the point with least 
	Face * f = total_min_face;
	// copy bfs code into this
	Array <Face *> array(500);
	queue<Face*> q;

	printf("Total_min = %lf\ntotal_median=%lftotal_min_face=%ld\n",total_min,total_median,(long)total_min_face);

	Face *temp;
	Face * t;

	while(!q.empty())
		q.pop();
	Bag<Face*> *bf = m->getFaces(); 
	Iterator<Face*> *iter = bf->StartIteration();
	while(Face* f2 = iter->GetNext())
	{
		f2->unMark();
	}
	bf->EndIteration(iter);

	q.push(f);
	int added_count =0;
	while(!q.empty())
	{
		temp = q.front();
		q.pop();

		// printf("%d\n",index);
		Edge *e = temp->getEdge();
		//printf("hi\n");
		//	   printf("%d\n",e);
		//Edge * temp1 = NULL;
		for (int counter =0; counter < 3; counter ++)
		{

			//  printf("hi2\n");
			//Edge *e1 = e->getOpposite();
			// printf("hi-inter\n");
			if(e->getOpposite()!= NULL)
			{
				// printf("hi3\n");
				//fflush(stdin);
				t = e->getOpposite()->getFace();
				// printf("hi4\n");
				if(!t->getMark() && hashcurvnew[(long)t] < total_median+0.01)
				{
					//  if(t==NULL)
					//  	   printf(":OOOOOOOOOOOOO\n");

					q.push(t);
					//	qindex.push(index-1);
					//					array.Add(t);
					t->setMark();
					added_count ++;

				}
			}
			//printf("test\n");
			e=e->getNext();
			//printf("test1\n");
		}
	}

	printf("I came till setting colors here %d\n",added_count);
	//for(;;);
	iter = bf->StartIteration();
	while(Face * f3 = iter->GetNext())
	{
		if(f3->getMark())
		{
			//f3->setColor(Vec3f(0.0,0.0,0.0));
			m->removeFace(f3);
		}
	}
	bf->EndIteration(iter);

	// for(int counter =0; counter < array.Count(); counter < 0)
	//	 {
	//	  		 //m->removeFace(array[counter]);
	//	  		 array[counter]->setColor(Vec3f(0.0,0.0,0.0));
	//	 }

}

int less_than(Vertex *a , Vertex * b)
{
	return (unsigned long)a < (unsigned long)b;
}

double getAngle(Vertex *a, Vertex*mid, Vertex *b)
{
	Vec3f p = a->get()-mid->get();
	Vec3f q = b->get()-mid->get();
	if(acos(p.Dot3(q)/p.Length()/q.Length())-M_PI>1e-2)
	{
		printf("acos is greater than pi?! %lf\n",acos(p.Dot3(q)/p.Length()/q.Length()));
		int unused = scanf("%*d");
    unused++;
	}

	float answer = 180.0/M_PI*acos(p.Dot3(q)/p.Length()/q.Length());
	if(answer -180>1e-2)
	{
		printf("did acos go above pi?\n");
		printf("acos is greater than pi?! %lf\n",acos(p.Dot3(q)/p.Length()/q.Length()));
		printf("my M_PI is %lf\n",M_PI);
		printf("my answer is %f %lf %lf %lf",answer,180.0,M_PI,acos(p.Dot3(q)/p.Length()/q.Length()));
		int unused = scanf("%*d");
    unused++;
	}
	return answer;
}

bool checkReflex(Vertex *a, Vertex*b, Vertex*c, Vertex*d)
{
	Vec3f p = c->get()-d->get();
	Vec3f q = a->get()-d->get();
	Vec3f res1;
	Vec3f::Cross3(res1,p,q);
	p = c->get()-d->get();
	q = b->get()-d->get();
	Vec3f res2;
	Vec3f::Cross3(res2,p,q);
	if((res1.Dot3(res2)/res1.Length()/res2.Length())<epsilonarea)
		return false;
	return true;
}


void findAllAngles(Face * f, float arr[3])
{
	double amin=10000000, amax=-100000;
	double aavg=0;
	int avg_count =0;
	Edge *e1 = f->getEdge();
	for(int count = 0 ; count < 3; count ++)
	{

		Edge * e2 = e1->getNext();
		Vec3f v1 = (*e1)[1]->get()-(*e1)[0]->get();
		Vec3f v2 = (*e2)[0]->get()-(*e2)[1]->get();
		float angle1 = 180.0/M_PI*(acos(v1.Dot3(v2)/v1.Length()/v2.Length()));

		if(amin > angle1)
			amin = angle1;
		if(amax < angle1)
			amax = angle1;
		aavg = aavg + angle1;
		avg_count ++;
		e1 = e1->getNext();
	}
	arr[0]=amin;
	arr[1]=aavg/avg_count;
	arr[2]=amax;
	//	printf("Individual Angle limits Amin : %lf   Aavg : %lf   Amax : %lf\n",amin,aavg/avg_count,amax);
}
void findAllAngles()
{
	Bag<Face *>* bf = m->getFaces();
	Iterator<Face *> * iter = bf->StartIteration();
	double amin=10000000, amax=-100000,aavg=0;
	int avg_count = 0;
	int cset = 0;
	while(Face * f = iter->GetNext())
	{
		Edge *e1 = f->getEdge();
		bool skewed = false;
		for(int count = 0 ; count < 3; count ++)
		{

			Edge * e2 = e1->getNext();
			Vec3f v1 = (*e1)[1]->get()-(*e1)[0]->get();
			Vec3f v2 = (*e2)[0]->get()-(*e2)[1]->get();
			float angle1 = 180.0/M_PI*(acos(v1.Dot3(v2)/v1.Length()/v2.Length()));

			if(angle1 > 120)
				skewed = true;
			if(amin > angle1)
				amin = angle1;
			if(amax < angle1)
				amax = angle1;
			aavg = aavg + angle1;
			avg_count ++;
			e1 = e1->getNext();
		}
		if(skewed)
		{
			cset ++;
			f->setColor(Vec3f(1.0,0,0));
		}
	}
	aavg = aavg/avg_count;
	bf->EndIteration(iter);
	// find aavg just to verify that it comes to 60
	printf("Angle limits Amin : %lf   Aavg : %lf   Amax : %lf cset = %d\n",amin,aavg,amax,cset);

}

double getAngle(Edge *a, Edge *b)
{
	Vec3f v1 = (*a)[1]->get()-(*a)[0]->get();
	Vec3f v2 = (*b)[1]->get()-(*b)[0]->get();

	return 180/M_PI*acos(v1.Dot3(v2)/v1.Length()/v2.Length());

}
void edgeSwap()
{
	Bag<Edge *> * be = m->getEdges();
	Iterator<Edge *> *iter = be->StartIteration();
	int random_test=0;
	while( Edge * e = iter->GetNext())
	{
		random_test++;
		//	if(random_test >50000)
		//	break;
		if(random_test %5000 ==0)
		{
			printf("%d\r",random_test);
			//	sleep(0.0001);
		}

		if(e->getOpposite()==NULL)
		{
			//	printf("I'm struck here \n");
			continue;
			for(;;);
			be->EndIteration(iter);
			return;
		}

		float arr1[3],arr2[3];
		findAllAngles(e->getFace(),arr1);
		findAllAngles(e->getOpposite()->getFace(),arr2);
		if(arr1[2] > 120 || arr2[2] >120)
		{
			if(getAngle(e,e->getNext())<100 || getAngle(e,e->getPrev()) <100 || getAngle(e->getOpposite(),e->getOpposite()->getNext())< 100 ||getAngle(e->getOpposite(),e->getOpposite()->getPrev())<100)
			{
				continue;
			}
			Vec3f vn1 = e->getFace()->getNormal();
			Vec3f vn2 = e->getOpposite()->getFace()->getNormal();
			if(vn2.Length()<epsilon|| vn1.Length()<epsilon)
			{
				continue;
				//be->EndIteration(iter);
				//return;
			}
			if( acos(vn1.Dot3(vn2)/vn1.Length()/vn2.Length())<15*M_PI/180.0)
			{
				Face * f = e->getFace();
				Face * fop = e->getOpposite()->getFace();
				Vertex * v[4];
				v[0]= (*e)[0];
				v[1]= (*e)[1];
				v[2]= (*(e->getNext()))[1];
				v[3]= (*(e->getOpposite()->getNext()))[1];
				//	printf("I did %d\n",rand());
				m->removeFace(f);
				m->removeFace(fop);
				m->addFace(v[3],v[1],v[2],Vec3f(1,1,1),Vec3f(1.0,0.0,0.0));//color3
				m->addFace(v[3],v[2],v[0],Vec3f(1,1,1),Vec3f(1.0,1.0,0.0));
			}
			else
			{
				e->getFace()->setEmit(Vec3f(1.0,1.0,1.0));
				e->getOpposite()->getFace()->setEmit(Vec3f(0.0,0.0,0.0));
			}
		}
	}
	be->EndIteration(iter);
}


bool collapse_old(Edge * e,int n)
{
	Vec3f	trycolors[10];
	trycolors[0].Set(0,1,1);
	trycolors[1].Set(1,0,1.0);
	trycolors[2].Set(1,1.0,0);
	trycolors[3].Set(1,0,0);
	trycolors[4].Set(0,1,0);
	trycolors[5].Set(0,0,1);




	Face *f1=NULL,*f2=NULL;
	f1 = e->getFace();
	if( e->getOpposite()!=NULL)
		f2 = e->getOpposite()->getFace();
	if(f2==NULL)
	{
		printf("because of this??\n");
		return false;
	}
	vector<Vertex *> vedge;
	vector<Vertex *> vedge1;
	vedge1.clear();
	vedge.clear();
	int bot = -1;
	Edge *temp = e->getNext()->getNext();
	Vertex *vbase = e->getVertex();
	Vertex *vnext = e->getNext()->getVertex();
	assert(vnext);
	assert(vbase);
	Vec3f mid  = (vnext->get()+vbase->get())/2.0;
	Vec3f pos = vnext->get();
	int pick = 116992;
	//float temparea = f1->getArea();
	Edge * store;
	while(1)
	{

		vedge1.push_back(temp->getVertex());
		if(n==pick && vedge1.size()==2)
			store = temp;
		if(temp->getOpposite()==NULL)
			return false;
		temp = temp->getOpposite()->getNext()->getNext();
		if(temp->getVertex()==vnext)
			break;

	}
	//	if(n==pick)
	//		m->removeFace(store->getFace());
	bot = vedge1.size()-1;
	temp = temp->getNext()->getNext()->getOpposite()->getNext()->getNext();
	while(1)
	{
		vedge1.push_back(temp->getVertex());
		temp = temp->getOpposite()->getNext()->getNext();
		if(temp->getOpposite()==NULL)
		{
			printf("because of this one? %d\n", rand());
			return false;
		}
		if(temp->getOpposite()->getNext()->getNext()->getVertex()==vbase)
		{
			break;
		}
	}
	for(int counterx = 0; counterx < bot; counterx++)
	{
		if(n==pick)
		{
			printf("angle was %lf\n",getAngle(vnext,vedge1[counterx],vedge1[counterx+1]));
		}
		if(getAngle(vnext,vedge1[counterx],vedge1[counterx+1])>150 || getAngle(vnext,vedge1[counterx+1],vedge1[counterx])>150)
		{
			returned_count++;
			printf("I returned due to obtuse angle constraint %lf %lf\t",getAngle(vnext,vedge1[counterx],vedge1[counterx+1]),getAngle(vnext,vedge1[counterx+1],vedge1[counterx]));
			return false;			
		}
		if(!checkReflex(vnext,vbase,vedge1[counterx],vedge1[counterx+1]))
		{
			returned_count++;
			printf("I returned due to reflex angle constraint\t");
			return false;
		}
	}

	sort(vedge1.begin(),vedge1.end(),less_than);
	vector<Vertex*>::iterator end1 = vedge1.end();
	vector<Vertex*>::iterator new_end = unique(vedge1.begin(),vedge1.end());
	if(new_end!= end1)
	{
		int d = end1-new_end;
		printf("Difference was %d \n",d);
		returned_count ++;

		return false;
	}
	if(n==pick)
	{
		printf("bot == %d\n",bot);
	}
	non_returned_count++;

	temp = e->getNext()->getNext();
	vedge.clear();

	while(1)
	{

		vedge.push_back(temp->getVertex());
		//temp->getVertex()->get().Write();
		Face * fdel = temp->getFace();

		if(temp->getOpposite())
			temp = temp->getOpposite()->getNext()->getNext();
		else
			return false;
		if(n!=pick)
			m->removeFace(fdel);
		else
			fdel->setColor(trycolors[vedge.size()]);

		if(temp->getVertex()==vnext)
		{
			if(n!=pick)
				m->removeFace(temp->getFace());
			else
				fdel->setColor(trycolors[vedge.size()]);
			break;
		}


		//	DEBUG("CAUGHT IN THE LOOP\n");
	}


	for(unsigned int counter =0; counter < vedge.size()-1; counter ++)
	{
		if(n==pick)
		{
			//continue;
			printf("%ld %ld %ld\n",(long int)vedge[counter], (long int)vedge[counter+1], (long int)vnext);
			vedge[counter]->get().Write();
			vedge[counter+1]->get().Write();
			vnext->get().Write();
			if(counter >1)
			{
				//	m->addFace(vedge[counter],vedge[counter+1],vnext,trycolors[counter],Vec3f(0.0,0,1));
			}
			else if(counter ==0 )
			{
				//		vedge[counter]->set(vedge[counter]->get()+Vec3f(0.5,0.0,0.0));
				//	m->addFace(vedge[counter],vedge[counter+1],vnext,Vec3f(1,0,1),Vec3f(0.0,0,1));
				//	printf("This the area %lf f1 was %lf\n",f23->getArea(),temparea);
			}
		}
		//	else if(n>200)
		//		m->addFace(vedge[counter],vedge[counter+1],vnext,Vec3f(1,0,0),Vec3f(0.0,0,1));
		else
			m->addFace(vedge[counter],vedge[counter+1],vnext,Vec3f(0,1,0),Vec3f(0.0,0,1));

	}
	//	m->getVertices()->Remove(vbase);
	vnext->set(mid);
	return true;
}
bool collapse(Edge * e,int n)
{
	Vec3f	trycolors[10];
	trycolors[0].Set(0,1,1);
	trycolors[1].Set(1,0,1.0);
	trycolors[2].Set(1,1.0,0);
	trycolors[3].Set(1,0,0);
	trycolors[4].Set(0,1,0);
	trycolors[5].Set(0,0,1);




	Face *f1=NULL,*f2=NULL;

	f1 = e->getFace();
	if( e->getOpposite()!=NULL)
		f2 = e->getOpposite()->getFace();
	if(f2==NULL)
	{
		//	printf("because of this??\n");
		return false;
	}
	vector<Vertex *> vedge;
	vector<Vertex *> vedge1;
	vedge1.clear();
	vedge.clear();
	int bot = -1;
	Edge *temp = e->getPrev();
	Vertex *vbase = (*e)[0];
	Vertex *vnext = (*e)[1];
	assert(vnext);
	assert(vbase);
	Vec3f mid  = (vnext->get()+vbase->get())/2.0;
	Vec3f pos = vnext->get();
	int pick = 116000992;
	//float temparea = f1->getArea();
	Edge * store;
	Vec3f nor=Vec3f(0,0,0);
	Vec3f t;

	double K=0;
	Face * f123;
	Vec3f vran0,vran1,vran2;
	while(1)
	{
		//assert((*temp)[1]==vbase);
		nor = nor - temp->getFace()->getNormal();
		Face * f123= temp->getFace();
		Vec3f vran0 = (*f123)[0]->get();
		Vec3f vran1 = (*f123)[1]->get();
		Vec3f vran2 = (*f123)[2]->get();
		K = K + det3x3s(vran0.x(),vran0.y(),vran0.z(),vran1.x(),vran1.y(),vran1.z(),vran2.x(),vran2.y(),vran2.z());

		vedge1.push_back(temp->getVertex());
		if(n==pick && vedge1.size()==2)
			store = temp;
		if(temp->getOpposite()==NULL)
			return false;
		temp = temp->getOpposite()->getPrev();
		if(temp->getVertex()==vnext)
		{
			nor = nor - temp->getFace()->getNormal();

			f123= temp->getFace();
			Vec3f vran0 = (*f123)[0]->get();
			Vec3f vran1 = (*f123)[1]->get();
			Vec3f vran2 = (*f123)[2]->get();
			K = K + det3x3s(vran0.x(),vran0.y(),vran0.z(),vran1.x(),vran1.y(),vran1.z(),vran2.x(),vran2.y(),vran2.z());

			break;
		}

	}
	//	printf("Finished loading the first part\n");
	//	fflush(stdout);

	//	if(n==pick)
	//		m->removeFace(store->getFace());
	bot = vedge1.size()-1;
	if(temp->getPrev()->getOpposite()==NULL)
		return false;

	temp = temp->getPrev()->getOpposite()->getPrev();
	while(1)
	{
		//assert((*temp)[1]==vnext);

		nor = nor - temp->getFace()->getNormal();
		f123= temp->getFace();
		Vec3f vran0 = (*f123)[0]->get();
		Vec3f vran1 = (*f123)[1]->get();
		Vec3f vran2 = (*f123)[2]->get();
		K = K + det3x3s(vran0.x(),vran0.y(),vran0.z(),vran1.x(),vran1.y(),vran1.z(),vran2.x(),vran2.y(),vran2.z());

		vedge1.push_back(temp->getVertex());
		if(temp->getOpposite()==NULL)
			return false;
		temp = temp->getOpposite()->getPrev();
		if(temp->getOpposite()==NULL)
			return false;
		if(temp->getOpposite()->getPrev()->getVertex()==vbase)
		{
			nor = nor - temp->getFace()->getNormal();
			f123= temp->getFace();
			Vec3f vran0 = (*f123)[0]->get();
			Vec3f vran1 = (*f123)[1]->get();
			Vec3f vran2 = (*f123)[2]->get();
			K = K + det3x3s(vran0.x(),vran0.y(),vran0.z(),vran1.x(),vran1.y(),vran1.z(),vran2.x(),vran2.y(),vran2.z());

			break;
		}
	}

	Vec3f vfinalpos=vnext->get();
	if(nor.Length()>1e-2)
	{
		K = K/nor.Length();
		nor.Normalize();
		double lambda = nor.Dot3(mid)-K;
		vfinalpos = mid - lambda*nor;
	}
	//	for(int counter =0; counter <=bot; counter++)
	//	{
	//		
	//	}
	//	printf("Finished loading the second part\n");
	//	fflush(stdout);
	for (int counterx = 0; counterx < bot; counterx++)
	{
		if(n==pick)
		{
			printf("angle was %lf\n",getAngle(vnext,vedge1[counterx],vedge1[counterx+1]));
		}
		if(getAngle(vnext,vedge1[counterx],vedge1[counterx+1])>150 || getAngle(vnext,vedge1[counterx+1],vedge1[counterx])>150)
		{
			returned_count++;
			//printf("I returned due to obtuse angle constraint\t");
			return false;			
		}
		if(!checkReflex(vnext,vbase,vedge1[counterx],vedge1[counterx+1]))
		{
			returned_count++;
			//printf("I returned due to reflex angle constraint\t");
			return false;
		}
	}

	sort(vedge1.begin(),vedge1.end(),less_than);
	vector<Vertex*>::iterator end1 = vedge1.end();
	vector<Vertex*>::iterator new_end = unique(vedge1.begin(),vedge1.end());
	if(new_end!= end1)
	{
		//int d = end1-new_end;
		//		printf("Difference was %d \n",d);
		returned_count ++;
		return false;
	}
	if(n==pick)
	{
		printf("bot == %d\n",bot);
	}
	non_returned_count++;

	temp = e->getPrev();
	vedge.clear();
	//	printf("About to start deleting\n");
	//	fflush(stdout);
	while(1)
	{

		//	assert((*temp)[1]==vbase);
		vedge.push_back(temp->getVertex());
		//temp->getVertex()->get().Write();
		Face * fdel = temp->getFace();

		if(temp->getOpposite())
			temp = temp->getOpposite()->getPrev();
		else
			return false;
		//	if(n!=pick)
		m->removeFace(fdel);
		//	else
		//	fdel->setColor(trycolors[vedge.size()]);

		if(temp->getVertex()==vnext)
		{
			//	if(n!=pick)
			m->removeFace(temp->getFace());
			//	else
			//			fdel->setColor(trycolors[vedge.size()]);
			break;
		}


		//	DEBUG("CAUGHT IN THE LOOP\n");
	}


	for(unsigned int counter =0; counter < vedge.size()-1; counter ++)
	{
		if(n==pick)
		{
			//continue;
			//printf("%d %d %d\n",vedge[counter],vedge[counter+1],vnext);
			vedge[counter]->get().Write();
			vedge[counter+1]->get().Write();
			vnext->get().Write();
			if(counter >1)
			{
				//	m->addFace(vedge[counter],vedge[counter+1],vnext,trycolors[counter],Vec3f(0.0,0,1));
			}
			else if(counter ==0 )
			{
				//		vedge[counter]->set(vedge[counter]->get()+Vec3f(0.5,0.0,0.0));
				//	m->addFace(vedge[counter],vedge[counter+1],vnext,Vec3f(1,0,1),Vec3f(0.0,0,1));
				//	printf("This the area %lf f1 was %lf\n",f23->getArea(),temparea);
			}
		}
		//	else if(n>200)
		//		m->addFace(vedge[counter],vedge[counter+1],vnext,Vec3f(1,0,0),Vec3f(0.0,0,1));
		else
			m->addFace(vedge[counter],vedge[counter+1],vnext,Vec3f(0,1,0),Vec3f(0.0,0,1));

	}
	//	m->getVertices()->Remove(vbase);
	//vnext->set(mid);
	vnext->set(vfinalpos);
	return true;
}
void Decimate(int num)
{
	int count =0;
	Bag <Edge *> *be = m->getEdges();
	printf("num edges = %d\n",m->numEdges());
	//	 scanf("%*d");
	// printf("Came here\n");
	int continue_count  =0;
	Iterator<Edge *> *iter = be->StartIteration();
	//int new_count=0;
	double lenlimit = 0.1;
	int confidence_continue_count = 0;
	double confidence_limit = 5;
	while(count < num)
	{
		// pick a random half edge and try to collapse
		Edge * random_edge = iter->GetNext();//be->ChooseRandom();
		//  sleep(0.01);

		if(random_edge==NULL)
		{
			//	printf("here\n");
			//	for(;;);
			be->EndIteration(iter);
			//edgeSwap();
			iter = be->StartIteration();
			//	printf("I ended up here\n");
			continue;
		}
		//	printf("I ended up here actually\n");
		//		  if(hash_decimate.find((int)random_edge)==hash_decimate.end())
		//		  {
		//		  		  hash_decimate[(int)random_edge]=1;
		//		  		  new_count++;
		//		  		  printf("%d %d %d\n",count,new_count,m->numEdges());
		//		  }
		//		  else
		{
			//	printf("%d %d %d\n",count,new_count,m->numEdges());
		}
		//	  else
		//	  {
		//			continue_count ++;
		//			printf("continuing\n");
		//			if(continue_count > 2000)
		//				break;
		//	  		  continue;
		//		}
		Vec3f rand_pos = random_edge->getVertex()->get();
		if(find_max_around_current_point(rand_pos)<confidence_limit)
		{
			confidence_continue_count++;
			if(confidence_continue_count > 0.1*m->numEdges())
			{
				//	printf("\nI did change it once: confidence limit is %lf\n",confidence_limit);
				confidence_continue_count=0;
				confidence_limit--;
			}
			continue;
		}

		//		double lmin,lmedian,lmax;
		//		double lvalue = find_max_around_current_point(rand_pos);
		//		lvalue = ((lvalue > lmax)?lmax:((lvalue<lmin)?lmin:lvalue));
		//		lvalue = (lvalue - lmin)/(lmax-lmin);
		//		lvalue = pow(lvalue,1);
		//		if((rand()%100000)/100000.0>lvalue)
		//		{
		//			continue;
		//		}
		if(getLength(random_edge)>lenlimit)
		{
			continue_count ++;
			if(continue_count%1000==0)
				//printf("%d\r",continue_count);
				if(continue_count > 2* m->numEdges())
				{
					//printf("\ncurrent lenlimit is %lf\n",lenlimit+0.1);
					lenlimit+=0.2;
					continue_count = 0;
				}
				continue;
		}
		//
		//	if(count <num/5)
		//	{
		//		if(getLength(random_edge)>0.5)
		//		{
		//			continue;
		//		}
		//	}
		//	else if( count < num/3)
		//	{
		//		if(getLength(random_edge)>1.0)
		//		{
		//			continue;
		//		}
		//	}
		//	else if(count < num/2)
		//	{
		//		if(getLength(random_edge)>1.5)
		//		{
		//			continue;
		//		}
		//	}



		count ++;
		// if(count >=7 && count<=7)
		//	  {
		//for(;;);
		//	printf("hi start\n");
		//sleep(2);
		//			printf("Calling collape_edge\n");
		bool success= collapse(random_edge,count);
		//		 printf("Returned from collapse_edge\n");
		if(!success)
		{
			//printf("I keep coming here %d\r",rand());
			count--;
		}
		//if(count%50==0)
		//	printf("%d\r",count);

		//printf("count %d\n",count);
		//sleep(0.1);
	}
	//sleep(2);

	// }
	be->EndIteration(iter);
}

//
//double getAngle(Edge *a, Edge *b)
//{
//	Vec3f v1 = (*a)[1]->get()-(*a)[0]->get();
//	Vec3f v2 = (*b)[1]->get()-(*b)[0]->get();
//	
//	return 180/M_PI*acos(v1.Dot3(v2)/v1.Length()/v2.Length());
//	
//}
//void edgeSwap()
//{
//	Bag<Edge *> * be = m->getEdges();
//	Iterator<Edge *> *iter = be->StartIteration();
//	int random_test=0;
//	while( Edge * e = iter->GetNext())
//	{
//		random_test++;
//	//	if(random_test >50000)
//	//	break;
//		if(random_test %50 ==0)
//		{
//			printf("%d\r",random_test);
//			//	sleep(0.0001);
//		}
//		
//		if(e->getOpposite()==NULL)
//		{
//			printf("I'm struck here \n");
//			continue;
//			for(;;);
//			be->EndIteration(iter);
//			return;
//		}
//		
//		float arr1[3],arr2[3];
//		findAllAngles(e->getFace(),arr1);
//		findAllAngles(e->getOpposite()->getFace(),arr2);
//		if(arr1[2] > 120 || arr2[2] >120)
//		{
//			if(getAngle(e,e->getNext())<100 || getAngle(e,e->getPrev()) <100 || getAngle(e->getOpposite(),e->getOpposite()->getNext())< 100 ||getAngle(e->getOpposite(),e->getOpposite()->getPrev())<100)
//			{
//				continue;
//			}
//			Vec3f vn1 = e->getFace()->getNormal();
//			Vec3f vn2 = e->getOpposite()->getFace()->getNormal();
//			if(vn2.Length()<epsilon|| vn1.Length()<epsilon)
//			{
//					continue;
//					//be->EndIteration(iter);
//					//return;
//			}
//			if( acos(vn1.Dot3(vn2)/vn1.Length()/vn2.Length())<15*M_PI/180.0)
//			{
//				Face * f = e->getFace();
//				Face * fop = e->getOpposite()->getFace();
//				Vertex * v[4];
//				v[0]= (*e)[0];
//				v[1]= (*e)[1];
//				v[2]= (*(e->getNext()))[1];
//				v[3]= (*(e->getOpposite()->getNext()))[1];
//			//	printf("I did %d\n",rand());
//				m->removeFace(f);
//				m->removeFace(fop);
//				m->addFace(v[3],v[1],v[2],Vec3f(1,1,1),Vec3f(1.0,0.0,0.0));//color3
//				m->addFace(v[3],v[2],v[0],Vec3f(1,1,1),Vec3f(1.0,1.0,0.0));
//			}
//			else
//			{
//				e->getFace()->setEmit(Vec3f(1.0,1.0,1.0));
//				e->getOpposite()->getFace()->setEmit(Vec3f(0.0,0.0,0.0));
//			}
//		}
//	}
//	be->EndIteration(iter);
//}

void checkWellFormed()
{
	//	return;
	printf("Entering checkWellFormed\n");
	Bag<Face *>*bf = m->getFaces();
	Iterator<Face*>*iter = bf->StartIteration();
	Array<Face*> arr(100);
	while(Face * f = iter->GetNext())
	{
		bool flag = false;
		arr.Clear();
		//BFS(f,1,&arr);
		Edge *e = f->getEdge();
		for(int counter =0; counter < 3; counter++)
		{
			if(e->getOpposite()==NULL)
			{
				continue;
				//printf("there was no opposite edge\n");
			}
			arr.Add(e->getOpposite()->getFace());
			e=e->getNext();
		}
		Vec3f vn = Vec3f(0.0,0.0,0.0) - f->getNormal();
		//printf("Entering Loop\n");
		fflush(stdout);
		for(int counter = 0; counter < arr.Count(); counter ++)
		{
			Vec3f vn1 = Vec3f(0,0,0) - arr[counter]->getNormal();
			double angle = 180.0/M_PI*acos(vn1.Dot3(vn)/vn1.Length()/vn.Length());
			if(angle > 150 )
			{
				printf("I found a point where angle = %lf  > 100 \n",angle);
				flag = true;
			}
		}
		//	printf("Leaving Loop\n");
		fflush(stdout);


		if(flag)
			f->setColor(Vec3f(1.0,0.0,0.0));
		else
			f->setColor(Vec3f(0,0,1));
		f->setColor(Vec3f(0,1,0));
		f->setEmit(Vec3f(0,1,0));

	}
	printf("Leaving checkWellFormed\n");
	bf->EndIteration(iter);
}

bool vlessthan(Vertex *a,Vertex*b)
{
	return a->getIndex()<b->getIndex();
}

void writePLY( char * filename )
{
	Array<Vertex*>* av = m->getVertices();
	vector <Vertex*> vv;
	Bag<Face*> * bf = m->getFaces();
	Iterator<Face *>* iter = bf->StartIteration();

	for(int counter =0; counter < av->Count(); counter++)
	{
		vv.push_back((*av)[counter]);
	}
	sort(vv.begin(),vv.end(),vlessthan);

	FILE*fp = fopen(filename,"w");
	if(fp==NULL)
	{
		printf("Couldn't open %s for writing\n",filename);
		bf->EndIteration(iter);
		return;
	}
	fprintf(fp,"ply\n");
	fprintf(fp,"format ascii 1.0\n");
	fprintf(fp,"comment test comment by author\n");
	fprintf(fp,"element vertex %d\n",(int)vv.size());
	fprintf(fp,"property float32 x\n");
	fprintf(fp,"property float32 y\n");
	fprintf(fp,"property float32 z\n");
	fprintf(fp,"element face %d\n",bf->Count());
	fprintf(fp,"property list uint8 float32 vertex_index\n");
	fprintf(fp,"end_header\n");
	for(unsigned int counter =0; counter < vv.size(); counter++)
	{
		fprintf(fp,"%lf %lf %lf\n",vv[counter]->x(),vv[counter]->y(),vv[counter]->z());
	}

	while(Face * f = iter->GetNext())
	{
		Edge *e1 = f->getEdge();
		Edge *etemp = e1;
		fprintf(fp,"3");
		do
		{
			int d = lower_bound(vv.begin(),vv.end(),(*etemp)[0],vlessthan)-vv.begin();

			fprintf(fp," %d",d);
			etemp=etemp->getNext();
		}while(e1!=etemp);
		fprintf(fp,"\n");
	}
	bf->EndIteration(iter);
	fclose(fp);
	printf("Finished writing to %s\n",filename);
}

void writeOBJ(const char * filename )
{
	Array<Vertex*>* av = m->getVertices();
	vector <Vertex*> vv;
	for(int counter =0; counter < av->Count(); counter++)
	{
		vv.push_back((*av)[counter]);
	}
	sort(vv.begin(),vv.end(),vlessthan);

	FILE*fp = fopen(filename,"w");
	if(fp==NULL)
	{
		printf("Couldn't open %s for writing\n",filename);
		return;
	}
	for(unsigned int counter =0; counter < vv.size(); counter++)
	{
		fprintf(fp,"v %lf %lf %lf\n",vv[counter]->x(),vv[counter]->y(),vv[counter]->z());
	}
	Bag<Face*> * bf = m->getFaces();
	Iterator<Face *>* iter = bf->StartIteration();
	while(Face * f = iter->GetNext())
	{
		Edge *e1 = f->getEdge();
		Edge *etemp = e1;
		fprintf(fp,"f");
		do
		{
			int d = lower_bound(vv.begin(),vv.end(),(*etemp)[0],vlessthan)-vv.begin();			
			fprintf(fp," %d",d+1);
			etemp=etemp->getNext();
		}while(e1!=etemp);
		fprintf(fp,"\n");
	}
	bf->EndIteration(iter);
	fclose(fp);
	printf("Finished writing to %s\n",filename);
}

void LoadOBJ(char * fname)
{
	FILE *fp = fopen(fname,"r");
	bool once = true;
	Array<Vertex*> *array;
	maxx = -1;
	maxy = -1;
	maxz = -1;
	while(!feof(fp))
	{
		double a,b,c;
		char ch;
		if( fscanf(fp,"%c",&ch) == EOF )
      {
      cerr << "EOF encountered in fscanf" << endl;
      }
		if(ch=='v')
		{
			if( fscanf(fp," %lf %lf %lf\n",&a,&b,&c) == EOF )
        {
        cerr << "EOF encountered in fscanf" << endl;
        }
			//	printf("vertex %lf %lf %lf\n",a,b,c);
			maxx = MAX(maxx,a);
			maxy = MAX(maxy,b);
			maxz = MAX(maxz,c);
			m->addVertex(Vec3f(a,b,c));
		}
		else if(ch=='f')
		{
			int p,q,r;
			if( fscanf(fp," %d %d %d\n",&p,&q,&r) == EOF )
        {
        cerr << "EOF encountered in fscanf" << endl;
        }
			//	printf("Face %d %d %d\n",p,q,r);
			if(once)
			{
				once = false;
				array = m->getVertices();
				m->addFace((*array)[p-1],(*array)[q-1],(*array)[r-1],Vec3f(1,1,1),color1);
			}
			else
			{
				m->addFace((*array)[p-1],(*array)[q-1],(*array)[r-1],Vec3f(1,1,1),color1);
			}
		}
		else
		{
			//			printf("praaalam\n");
			//			for(;;);
		}
	}
	fclose(fp);
}

void LoadOBJ(FILE *fp, Mesh *mesh)
{

	bool once = true;
	Array<Vertex*> *array;
	while(!feof(fp))
	{
		double a,b,c;
		char ch;
		if ( fscanf(fp,"%c",&ch) == EOF ) 
      {
      cerr << "EOF encountered in fscanf" << endl;
      }
		if(ch=='v')
		{
			if( fscanf(fp," %lf %lf %lf\n",&a,&b,&c) == EOF )
        {
        cerr << "EOF encountered in fscanf" << endl;
        }
			//	printf("vertex %lf %lf %lf\n",a,b,c);
			mesh->addVertex(Vec3f(a,b,c));
		}
		else if(ch=='f')
		{
			int p,q,r;
			if( fscanf(fp," %d %d %d\n",&p,&q,&r) == EOF )
        {
        cerr << "EOF encountered in fscanf" << endl;
        }
			//	printf("Face %d %d %d\n",p,q,r);
			if(once)
			{
				once = false;
				array = mesh->getVertices();
				mesh->addFace((*array)[p-1],(*array)[q-1],(*array)[r-1],Vec3f(1,1,1),color1);
			}
			else
			{
				mesh->addFace((*array)[p-1],(*array)[q-1],(*array)[r-1],Vec3f(1,1,1),color1);
			}
		}
		else
		{
			//			printf("praaalam\n");
			//			for(;;);
		}
	}
	fclose(fp);
	printf("Finished loading obj file\n");
}

int generate_hash(Face *f, int xbins = 0, int ybins = 0, int zbins = 0, double xmin = 0, double ymin = 0, double zmin = 0, double xmax =0, double ymax = 0, double zmax = 0)
{
	static int xb,yb,zb;
	static double minx,miny,minz,maxx,maxy,maxz;
	if(xbins >0)
	{
		xb = xbins;
		yb = ybins;
		zb = zbins;
		minx = xmin;
		miny = ymin;
		minz = zmin;
		maxx = xmax;
		maxy = ymax;
		maxz = zmax;
		return -1;
	}
	Vec3f center = f->getCenter();

	int num1 = (center.x()-minx)/(maxx-minx)*xb;
	int num2 = (center.y()-miny)/(maxy-miny)*yb;
	int num3 = (center.z()-minz)/(maxz-minz)*zb;

	// cheating! ;)
	if(xbins == -1)
		return num1;
	if(xbins == -2)
		return num2;
	if(xbins == -3)
		return num3;

	int total = num1+num2*xb+num3*xb*yb;
	return total;
}
double find_closest(Face *f, vtksys::hash_multimap<int,Face*> &hash,int xbins = 0, int ybins =0, int zbins =0)
{
	static int xb,yb,zb;
	if(xbins!=0)
	{
		xb = xbins;
		yb = ybins;
		zb = zbins;
		return -1;
	}
	int n[3];
	n[0] = generate_hash(f,-1);// the cheating trick
	n[1] = generate_hash(f,-2);
	n[2] = generate_hash(f,-3);



	int depth = 0;
	double mindist = 1000000;
	double dist;
	int number = 0;
	while(1)
	{
		//printf("current depth = %d %d %d %d\n",depth,xb,yb,zb);
		for(int cox = -depth; cox<=depth; cox++)
		{
			for(int coy = -depth; coy<=depth;coy++)
			{
				for(int coz = -depth; coz<=depth; coz++)
				{

					if(abs(cox)==depth || abs(coy)==depth || abs(coz)==depth)
					{
						number = (n[0]+cox)+(n[1]+coy)*xb+(n[2]+coz)*xb*yb;
						pair<vtksys::hash_multimap<int,Face*>::const_iterator,vtksys::hash_multimap<int,Face*>::const_iterator> pa = hash.equal_range(number);

						for(vtksys::hash_multimap<int,Face*>::const_iterator i = pa.first; i!=pa.second; ++i)
						{
							Face *f1 = (*i).second;
							dist = ((f1->getCenter())-(f->getCenter())).Length();
							if(dist < mindist)
								mindist = dist;
						}
					}
				}
			}
		}

		if(depth==0 || mindist > 100000)
		{
			depth = depth+1;
			continue;
		}
		/*	if(depth == 2)
		{
		printf("I did try to quit with depth == 2\n");
		}*/
		break;
	}
	return mindist;
}
void print_all_data(Mesh *mesh)
{
	printf("Num Faces : %d\n",mesh->numFaces());
	printf("Num Vertices : %d\n",mesh->numVertices());
	printf("Num Edges : %d\n",mesh->numEdges());

}
double find_distance_between_two_3D_files(char *fname1, char*fname2)
{
	int xbins = 50,ybins = 50, zbins = 20;
	FILE *fp1, *fp2;
	fp1 = fopen(fname1,"r");
	if(fp1==NULL)
	{
		printf("Error: File '%s' not found\n",fname1);
		return -1;
	}
	Mesh * m1, *m2;
	m1 = new Mesh();
	//loading obj made easy :) call and pray strategy. No checks for loading obj file.
	LoadOBJ(fp1,m1);
	fclose(fp1);
	fp2 = fopen(fname2,"r");
	if(fp2==NULL)
	{
		printf("Error: File '%s' not found\n",fname2);
		return -1;
	}

	m2 = new Mesh();
	LoadOBJ(fp2,m2);
	fclose(fp2);
	Vec3f min,max; 
	m1->getBoundingBox()->Get(min,max);
	print_all_data(m1);
	print_all_data(m2);
	// set up hashing function's static variables
	generate_hash(NULL,xbins,ybins,zbins,min.x(),min.y(),min.z(),max.x(),max.y(),max.z());

	// create a hashface for the two meshes

	vtksys::hash_multimap<int, Face*> hashface1,hashface2; // bow down to stl
	Bag<Face*>* bf = m1->getFaces();
	Iterator<Face*>*iter = bf->StartIteration();
	while(Face *f = iter->GetNext())
	{
		hashface1.insert(vtksys::hash_multimap<int,Face*>::value_type(generate_hash(f),f));
	}
	bf->EndIteration(iter);
	printf("Finished creating the hash for mesh 1\n");
	bf = m2->getFaces();
	iter = bf->StartIteration();
	double distance_sum = 0;
	find_closest(NULL,hashface1,xbins,ybins,zbins);
	int count = 0;
	printf("Starting the main loop:\n");
	while(Face *f = iter->GetNext())
	{
		//hashface2.insert(vtksys::hash_multimap<int,Face*>::value_type(generate_hash(f),f));
		distance_sum +=find_closest(f,hashface1);
		printf("%d %0.2lf\r",(int)distance_sum,count*1.0/m2->numFaces()*100.0);
		count++;
	}
	printf("\nTotal distance between the two models: %lf\nAverage Distance between the faces: %lf\n",distance_sum,distance_sum/m2->numFaces());
	bf->EndIteration(iter);
	//end of hash creation
	return distance_sum;
}

void calcnew(Vertex * a, char *buff)
{
	Vec3f verttemp = a->get();
	for(int counter =0; counter < 3; counter ++)
	{

		float d1 = verttemp[counter];
		int d = *((int*)&d1);
		for(int counter2 = 0 ; counter2 < 4; counter2 ++ )
		{
			buff[counter*4 + counter2] = d & 255;
			d = d >>8;
		}
	}
}




struct lface
{
	bool operator()(const Face* s1, const Face* s2) const
	{
		return s1 < s2;
	}
};

void scale_vertices(double a, double b, double c)
{
	current_scaling.Set(current_scaling.x()/a,current_scaling.y()/b,current_scaling.z()/c);
	Array<Vertex*> * arr = m->getVertices();
	Vec3f vt;
	for(int counter =0; counter < arr->Count(); counter++)
	{
		vt = (*arr)[counter]->get();
		vt.Set(vt[0]*a,vt[1]*b,vt[2]*c);
		(*arr)[counter]->set(vt);
	}
	return ; 
}

void delete_small_components(int num_to_delete = 40)
{
	printf("Entered delete_small_components\n");
	Bag<Face*> *bf = m->getFaces();
	Iterator<Face*> *iter = bf->StartIteration();
	Array<Face*> arr(200);
	int smallest=2343223;
	int num_deletes = 0;
	int num_faces =0;
	bool flag = true;
	while (Face*f = iter->GetNext())
	{
		arr.Clear();
		num_faces++;
		BFS(f,10,&arr);
		if(arr.Count() < smallest)
			smallest=arr.Count();
		if(num_faces%100)
			printf("%0.2lf%%\r",num_faces*1.0/m->numFaces()*100.0);
		if(arr.Count()<num_to_delete)
		{
			printf("I am going to delete %d faces\n",arr.Count());
			for(int counter=0; counter<arr.Count(); counter++)
			{
				//	if(arr[counter]!=f)
				m->removeFace(arr[counter]);
			}
			num_deletes++;
			if(num_faces*1.0/m->numFaces()>0.9 && flag )
			{
				flag = false;
				num_deletes = 0;
				num_faces=0;
				bf->EndIteration(iter);
				printf("Resetting now\n");
				iter = bf->StartIteration();
			}
		}
	}
	printf("Smallest = %d\n",smallest);
	bf->EndIteration(iter);
	printf("leaving delete_small_components\n");
}

//void make_hypothesis_model(void)
//{
//	int array1[][3]={ 0,0,0,
//					1,0,0,
//					0,1,0,
//					0,0,1,
//					1,0,1,
//					0,1,1,
//					1,1,0,
//					1,1,1};
//	int faces[][4]={0,1,4,3,
//					0,3,5,2,
//					3,4,7,5,
//					1,6,7,4,
//					6,2,5,7,
//					0,2,6,1
//	};
//	for(int counter=0; counter < 8; counter++)
//	{
//		m->addVertex(Vec3f(array1[counter][0],array1[counter][1],array1[counter][2]));
//	}
//	for(int counter=0; counter < 8; counter++)
//	{
//		m->addVertex(Vec3f(array1[counter][0]+5.0,array1[counter][1],array1[counter][2]));
//	}
//	Array<Vertex*>*array = m->getVertices();
//	printf("size %d\n",array->Count());
//	scanf("%*c");
//	for(int counter=0; counter<6;counter++)
//	{
//		m->addFace((*array)[faces[counter][0]],(*array)[faces[counter][1]],(*array)[faces[counter][2]],Vec3f(1.0,0.0,0.0),Vec3f(1.0,1.0,1.0));
//		m->addFace((*array)[faces[counter][0]],(*array)[faces[counter][2]],(*array)[faces[counter][3]],Vec3f(1.0,0.0,0.0),Vec3f(1.0,1.0,1.0));
//		m->addFace((*array)[faces[counter][0]+8],(*array)[faces[counter][1]+8],(*array)[faces[counter][2]+8],Vec3f(1.0,0.0,0.0),Vec3f(1.0,1.0,1.0));
//		m->addFace((*array)[faces[counter][0]+8],(*array)[faces[counter][2]+8],(*array)[faces[counter][3]+8],Vec3f(1.0,0.0,0.0),Vec3f(1.0,1.0,1.0));
//	}
//	m->addVertex(Vec3f(0,0,0.5));
//	m->addVertex(Vec3f(1,1,0.5));
//	m->addFace((*array)[],,Vec3f(1.0,0.0,0.0),Vec3f(1.0,1.0,1.0));
//}

void check_consistent()
{
	Bag<Face*> *bf = m->getFaces();
	Iterator<Face*> *iter = bf->StartIteration();
	while(Face *f = iter->GetNext())
	{
		Edge *e = f->getEdge();
		for(int counter =0; counter<3; counter++)
		{
			if(e==NULL||e->getOpposite()==NULL)
			{
				printf("its not proper");
				e->getFace()->setColor(Vec3f(1,1,1));
				e->getFace()->setEmit(Vec3f(1,1,1));
			}
			e=e->getNext();

		}
		f->setColor(Vec3f(1.0,1.0,1.0));
	}
	bf->EndIteration(iter);
	printf("Its proper!\n");
}
Vec3f get_color_code(int a)
{
	a = (a>int(lmax+0.5))?(int(lmax+0.5)):((a<int(lmin+0.5))?int(lmin+0.5):a);
	//printf("|%d|",a);
	return getColorCode(a,lmedian,lmin,lmax);
}
void find_min_median_max(void)
{
	vector<double> lvalues;
	lvalues.reserve(100000);
	lvalues.clear();
	for(int cx =0; cx <maxx; cx++)
		for(int cy =0; cy <maxy; cy++)
			for(int cz =0; cz <maxz; cz++)
			{
				if(pmatrix[cx][cy][cz]>0)
				{
					lvalues.push_back(pmatrix[cx][cy][cz]);
				}
			}
			sort(lvalues.begin(),lvalues.end());
			lmedian = lvalues[0.7*lvalues.size()];
			lmin = lvalues[0.05*lvalues.size()];
			lmax = lvalues[0.95*lvalues.size()];
			printf("I found %lf %lf %lf\n",lmin,lmedian,lmax);

}
int find_max_around_current_point(Vec3f &a)
{
	//	printf("current z scaling is %lf\n",current_scaling.z());
	int i=int(a.x()+0.5),j=int(a.y()+0.5),k=int(a.z()/current_scaling.z()+0.5);
	int max_curr = -1;
	for(int cx = i-2;cx<=i+2;cx++)
		for(int cy = j-2;cy<=j+2;cy++)
			for(int cz = k-1;cz<=k+1;cz++)
			{
				if(cx>=0 && cx <maxx && cy>=0 && cy<maxy && cz>=0 && cz<maxz)
					if(max_curr<pmatrix[cx][cy][cz])
					{
						max_curr=pmatrix[cx][cy][cz];
					}
			}
			return max_curr;
}
void Mesh::PaintConfidence()
{
	glDisable(GL_LIGHTING);
	float a,b;
	glGetFloatv(GL_POINT_SIZE_RANGE,&a);
	glGetFloatv(GL_POINT_SIZE_GRANULARITY,&b);
	// printf("Point size range allowed %lf , min and max being %lf %lf\n", a-b,b,a);

	glLineWidth(1);
	
	//
	//  glBegin (GL_POINTS);
	//	int counter =0;
	//  while (counter < vertices->Count()) {
	//		glColor3f(1,0,0);
	//		Vertex *v = (*vertices)[counter];
	//		counter ++;
	//    Vec3f a = v->get();
	//    Vec3f color = get_color_code(pmatrix[int(a.x()+0.5)][int(a.y()+0.5)][int(a.z()/2+0.5)]);
	//    glColor3f(color.r(),color.g(),color.b());
	//    glVertex3f(a.x(),a.y(),a.z());
	//  }
	//
	//  glEnd();
	glBegin(GL_TRIANGLES);
	Bag<Face *> *bf = m->getFaces();
	Iterator<Face*> *iter = bf->StartIteration();
	int L;
	Vec3f val,color;
	while(Face *f = iter->GetNext())
	{
		for(int counter=0; counter<3; counter++)
		{
			val = (*f)[counter]->get();
			L= find_max_around_current_point(val);
			//		if(L<pmatrix[int(val.x()+0.5)][int(val.y()+0.5)][int(val.z()/2+0.5)])
			//		printf("something is wrong\n");
			//			L = pmatrix[int(val.x()+0.5)][int(val.y()+0.5)][int(val.z()/2+0.5)];
			color = get_color_code(L);
			glColor3f(color.r(),color.g(),color.b());
			glVertex3f(val.x(),val.y(),val.z());
		}
	}
	glEnd();
	bf->EndIteration(iter);

	glEnable(GL_LIGHTING);

	HandleGLError(); 
}

void get_distance_statistics()
{
	Bag<Edge*>* be = m->getEdges();
	Iterator<Edge*>* iter = be->StartIteration();
	double mind = 1231232,maxd = -1,meand=0;
	int numd = 0;
	while(Edge * e = iter->GetNext())
	{
		numd++;
		double distance_e = ((*e)[0]->get()-(*e)[1]->get()).Length();
		if(mind >distance_e)
			mind = distance_e;
		if(maxd <distance_e)
			maxd = distance_e;
		meand+=distance_e;
		numd++;
	}
	printf("%lf %lf %lf\n",mind,meand/numd,maxd);
	be->EndIteration(iter);
}
int compare_decimations(char *filename1)
{
	int dec =0;
	int xbins = 50,ybins = 50, zbins = 20;

	Mesh *m1 = new Mesh();


	for(int x = 0;x <520;x++)
		for(int y = 0;y <520;y++)
			for(int z = 0;z <106;z++)
				pmatrix[x][y][z]=0; 



	m=m1;
	// "E:\\Arun\\My matlab codes\\deconvolved-lower3-higher8-window2";

	load_points_without_normal(filename1);

	marchtetra();

	check_consistent();

	scale_vertices(current_scaling.x(),current_scaling.y(),current_scaling.z());
	for(int counter =0; counter <3; counter++ )
	{
		if(0)
			SmoothSurface_complex(1);
		else
			SmoothSurface();
		printf("Initial smoothing %d time(s)\n", counter+1);
	}

	printf("Num vertices (before decimating) %d\n",m->numVertices());
	printf("Num Faces (before decimating) %d\n",m->numFaces());
	printf("dec?");
	//dec = 0;



	Mesh *m2 = new Mesh();


	for(int x = 0;x <520;x++)
		for(int y = 0;y <520;y++)
			for(int z = 0;z <106;z++)
				pmatrix[x][y][z]=0; 




	m=m2;
	// "E:\\Arun\\My matlab codes\\deconvolved-lower3-higher8-window2";


	load_points_without_normal(filename1);

	marchtetra();

	check_consistent();

	scale_vertices(current_scaling.x(),current_scaling.y(),current_scaling.z());
	for(int counter =0; counter <4; counter++ )
	{
		if(counter>2)
			SmoothSurface_complex(1);
		else
			SmoothSurface();
		printf("Initial smoothing %d time(s)\n", counter+1);
	}

	printf("Num vertices (before decimating) %d\n",m->numVertices());
	printf("Num Faces (before decimating) %d\n",m->numFaces());
	printf("dec?");
	dec = 0;
	//dec = 3*m->numFaces()/8;
	int temp_count =0;
	while(temp_count <10)
	{
		temp_count++;
		dec = m->numFaces()/4;
		Decimate(dec);
		printf("Num vertices (After decimating) %d\n",m->numVertices());
		printf("Num Faces (After decimating) %d\n",m->numFaces());
		SmoothSurface_complex(1);

		{
			Vec3f min,max; 
			m1->getBoundingBox()->Get(min,max);
			print_all_data(m1);
			print_all_data(m2);
			// set up hashing function's static variables
			generate_hash(NULL,xbins,ybins,zbins,min.x(),min.y(),min.z(),max.x(),max.y(),max.z());

			// create a hashface for the two meshes

			vtksys::hash_multimap<int, Face*> hashface1,hashface2; // bow down to stl
			Bag<Face*>* bf = m1->getFaces();
			Iterator<Face*>*iter = bf->StartIteration();
			while(Face *f = iter->GetNext())
			{
				hashface1.insert(vtksys::hash_multimap<int,Face*>::value_type(generate_hash(f),f));
			}
			bf->EndIteration(iter);
			printf("Finished creating the hash for mesh 1\n");
			bf = m2->getFaces();
			iter = bf->StartIteration();
			double distance_sum = 0;
			find_closest(NULL,hashface1,xbins,ybins,zbins);
			int count = 0;
			printf("Starting the main loop:\n");
			while(Face *f = iter->GetNext())
			{
				//hashface2.insert(vtksys::hash_multimap<int,Face*>::value_type(generate_hash(f),f));
				distance_sum +=find_closest(f,hashface1);
				if(((int)distance_sum)%100000==0 )
					printf("%0.2lf %0.2lf\n",distance_sum/count,count*1.0/m2->numFaces()*100.0);
				count++;
			}
			printf("\nTotal distance between the two models: %lf\nAverage Distance between the faces: %lf\n",distance_sum,distance_sum/m2->numFaces());
			bf->EndIteration(iter);
		}


		m = m2;
		m2 = m1;
		m1 = m;
		Vec3f min,max; 
		m1->getBoundingBox()->Get(min,max);
		print_all_data(m1);
		print_all_data(m2);
		// set up hashing function's static variables
		generate_hash(NULL,xbins,ybins,zbins,min.x(),min.y(),min.z(),max.x(),max.y(),max.z());

		// create a hashface for the two meshes
		{
			vtksys::hash_multimap<int, Face*> hashface1,hashface2; // bow down to stl
			Bag<Face*>* bf = m1->getFaces();
			Iterator<Face*>*iter = bf->StartIteration();
			while(Face *f = iter->GetNext())
			{
				hashface1.insert(vtksys::hash_multimap<int,Face*>::value_type(generate_hash(f),f));
			}
			bf->EndIteration(iter);
			printf("Finished creating the hash for mesh 1\n");
			bf = m2->getFaces();
			iter = bf->StartIteration();
			double distance_sum = 0;
			find_closest(NULL,hashface1,xbins,ybins,zbins);
			int count = 0;
			printf("Starting the main loop:\n");
			while(Face *f = iter->GetNext())
			{
				//hashface2.insert(vtksys::hash_multimap<int,Face*>::value_type(generate_hash(f),f));
				distance_sum +=find_closest(f,hashface1);
				if(((int)distance_sum)%100000==0 )
					printf("%0.2lf %0.2lf\n",distance_sum/count,count*1.0/m2->numFaces()*100.0);
				count++;
			}
			printf("\nTotal distance between the two models: %lf\nAverage Distance between the faces: %lf\n",distance_sum,distance_sum/m2->numFaces());
			bf->EndIteration(iter);
		}
		m = m2;
		m2 = m1;
		m1 = m;

		m = m2;
		char buff[1024];
		sprintf(buff,"meshing_decimated%d.obj",temp_count);
		writeOBJ(buff);
	}
	return 0;
}
int compare_npts_files(char *filename1, char*filename2)
{

	int dec =0;
	int xbins = 100,ybins = 100, zbins = 20;

	Mesh *m1 = new Mesh();


	for(int x = 0;x <520;x++)
		for(int y = 0;y <520;y++)
			for(int z = 0;z <106;z++)
				pmatrix[x][y][z]=0; 



	m=m1;
	// "E:\\Arun\\My matlab codes\\deconvolved-lower3-higher8-window2";

	load_points_without_normal(filename1);

	marchtetra();

	check_consistent();

	scale_vertices(current_scaling.x(),current_scaling.y(),current_scaling.z());
	for(int counter =0; counter <4; counter++ )
	{
		if(counter>2)
			SmoothSurface_complex(1);
		else
			SmoothSurface();
		printf("Initial smoothing %d time(s)\n", counter+1);
	}

	printf("Num vertices (before decimating) %d\n",m->numVertices());
	printf("Num Faces (before decimating) %d\n",m->numFaces());
	printf("dec?");
	dec = 0;
	//dec = 3*m->numFaces()/8;
	Decimate(dec);
	printf("Num vertices (After decimating) %d\n",m->numVertices());
	printf("Num Faces (After decimating) %d\n",m->numFaces());

	writeOBJ("mesh1.obj");
	Mesh *m2 = new Mesh();


	for(int x = 0;x <520;x++)
		for(int y = 0;y <520;y++)
			for(int z = 0;z <106;z++)
				pmatrix[x][y][z]=0; 




	m=m2;
	// "E:\\Arun\\My matlab codes\\deconvolved-lower3-higher8-window2";


	load_points_without_normal(filename2);

	marchtetra();

	check_consistent();

	scale_vertices(current_scaling.x(),current_scaling.y(),current_scaling.z());
	for(int counter =0; counter <4; counter++ )
	{
		if(counter>2)
			SmoothSurface_complex(1);
		else
			SmoothSurface();
		printf("Initial smoothing %d time(s)\n", counter+1);
	}

	printf("Num vertices (before decimating) %d\n",m->numVertices());
	printf("Num Faces (before decimating) %d\n",m->numFaces());
	printf("dec?");
	dec = 0;
	//dec = 3*m->numFaces()/8;
	Decimate(dec);
	printf("Num vertices (After decimating) %d\n",m->numVertices());
	printf("Num Faces (After decimating) %d\n",m->numFaces());
	writeOBJ("mesh2.obj");


	Vec3f min,max; 
	m1->getBoundingBox()->Get(min,max);
	print_all_data(m1);
	print_all_data(m2);
	// set up hashing function's static variables
	generate_hash(NULL,xbins,ybins,zbins,min.x(),min.y(),min.z(),max.x(),max.y(),max.z());

	// create a hashface for the two meshes

	vtksys::hash_multimap<int, Face*> hashface1,hashface2; // bow down to stl
	Bag<Face*>* bf = m1->getFaces();
	Iterator<Face*>*iter = bf->StartIteration();
	while(Face *f = iter->GetNext())
	{
		hashface1.insert(vtksys::hash_multimap<int,Face*>::value_type(generate_hash(f),f));
	}
	bf->EndIteration(iter);
	printf("Finished creating the hash for mesh 1\n");
	bf = m2->getFaces();
	iter = bf->StartIteration();
	double distance_sum = 0;
	find_closest(NULL,hashface1,xbins,ybins,zbins);
	int count = 0;
	printf("Starting the main loop:\n");
	while(Face *f = iter->GetNext())
	{
		//hashface2.insert(vtksys::hash_multimap<int,Face*>::value_type(generate_hash(f),f));
		distance_sum +=find_closest(f,hashface1);
		if(((int)distance_sum)%2000==0 )
			printf("%0.2lf %0.2lf\n",distance_sum/count,count*1.0/m2->numFaces()*100.0);
		count++;
	}
	printf("\nTotal distance between the two models: %lf\nAverage Distance between the faces: %lf\n",distance_sum,distance_sum/m2->numFaces());
	bf->EndIteration(iter);
	//end of hash creation
	return distance_sum;

}

#define DEBUGP 



Vec3f fit_circle(Array<Face*> *arr,float &fit,int max_iterations=10,double delta = 1.0)
{
	fit = INT_MAX/2;
	//DEBUGP("Entered fit_circle\n");
	if(arr->Count()==0)
	{
		return Vec3f(0,0,0);
	}
	//find the mean points and initialize it as the center thingy
	Vec3f mean = Vec3f(0,0,0);
	for(int counter=0; counter < arr->Count();counter++)
	{
		mean = mean + (*arr)[counter]->getCenter();
	}
	mean = mean/arr->Count();
	float r_est = 0;
	Vec3f c_est = mean;// begin with the mean as a good estimate;
	double sum_x,sum_y;

	//DEBUGP("Began simple gradient descent\n");
	//DEBUGP("Array->Count() = %d\n",arr->Count());
	// simple gradient descent 
	for (int iter =0;iter<max_iterations; iter++)
	{
		r_est =0;
		for(int co = 0; co<arr->Count();co++)
		{
			r_est = r_est + ((*arr)[co]->getCenter()-c_est).Length();
		}
		r_est = r_est/arr->Count();
		sum_x=0;
		sum_y=0;
		double di,weight;
		Vec3f new_est = Vec3f(0,0,0);
		//DEBUGP("Done finding r_est\n");
		for(int counter=0; counter<arr->Count(); counter++)
		{
			di = ((*arr)[counter]->getCenter()-c_est).Length();
			weight = (di-r_est)/di;
			new_est = new_est + (c_est-(*arr)[counter]->getCenter())*weight;
		}
		//if(new_est.Length()>2)
		{
			new_est.Normalize();
		}
		c_est = c_est - 2*delta*new_est;

		double J = 0;
		for(int counter=0; counter<arr->Count(); counter++)
		{
			J = J + pow((((*arr)[counter]->getCenter()-c_est).Length()-r_est),2);
		}
		J = J/arr->Count();
		if(J<0.005)
			break;
	}
	if(c_est.x()<0 || c_est.y()<0|| c_est.z()<0 || c_est.x()>=maxx || c_est.y()>=maxy || c_est.z()>=maxz)
		return Vec3f(-1,-1,-1);

	double J = 0;
	for(int counter=0; counter<arr->Count(); counter++)
	{
		J = J + pow((((*arr)[counter]->getCenter()-c_est).Length()-r_est),2);
	}
	J = J/arr->Count();
	if(J >5 || J!=J)
		return Vec3f(-1,-1,-1);
	fit = J;
	//printf("J=%lf\n",J);
	//DEBUGP("Returning from fit_circle\n");
	return c_est;
}

#define VOTES(a,b,c) votes[int(a)*maxy*maxz+int(b)*maxz+int(c)]
#define FIT(a,b,c) fit_values[int(a)*maxy*maxz+int(b)*maxz+int(c)]

float *fit_values;



//#define DISPFIT
void set_only_max_votes()
{
	bool flag;
	for(int cox =0; cox<maxx;cox++)
	{
		int bx = MAX(cox-2,0);
		int ex = MIN(cox+2,maxx-1);
		for(int coy=0;coy<maxy;coy++)
		{
			int by = MAX(coy-2,0);
			int ey = MIN(coy+2,maxy-1);
			for(int coz=0;coz<maxz;coz++)
			{
				int bz = MAX(coz-2,0);
				int ez = MIN(coz+2,maxz-1);
				flag = false;
				for(int cx = bx;cx<=ex;cx++)
					for(int cy=by;cy<=ey;cy++)
						for(int cz=bz;cz<ez;cz++)
						{
							if(!(cx==cox&&cy==coy&&cz==coz) && VOTES(cx,cy,cz)>VOTES(cox,coy,coz))
							{
								flag = true;
								goto label;
							}
						}
label:
						if(flag == true)
						{
							VOTES(cox,coy,coz)=0;
						}

			}
		}
	}
}

void get_centerlines_voxels();

void getSmoothed(vector<Vec3f> &line)
{
	int pc = 0;
	Vec3f temp;
	while(pc <2)
	{
		for(unsigned int counter=1; counter<line.size()-1; counter++)
		{	
			line[counter].Set((line[counter-1].x()+line[counter+1].x())/2.0,(line[counter-1].y()+line[counter+1].y())/2.0,(line[counter-1].z()+line[counter+1].z())/2.0);
		}
		pc++;
	}
}

void Mesh::PaintCenterLines()
{
	char filename_microgliatrace[] = "half_vessel_votes1_centerline.txt";

glLineWidth(5);
	{
		FILE *fp = fopen(filename_microgliatrace,"r");
		if(fp == NULL)
		{
			printf("Couldnt open %s file for reading\n",filename_microgliatrace);
			return ;
		}
		printf("Reading file %s\n",filename_microgliatrace);
		float x,y,z;
		bool start = true;
		vector<Vec3f> line;
		line.clear();

		while(fscanf(fp,"%f %f %f",&x,&y,&z)>0)
		{
			if(x==-1)
			{
			//	getSmoothed(line);

			//	printf("Starting to add trace\n");
			glDisable(GL_LIGHTING);
				glBegin(GL_LINE_STRIP);
				glColor3f(1,0,0);
				for(unsigned int counter=0; counter< line.size(); counter++)
				{
					glVertex3f(line[counter].x(),line[counter].y(),line[counter].z());
				}
				glEnd();
				glEnable(GL_LIGHTING);
				line.clear();
				start = true;
				continue;
			}
			//fscanf(fp,"%*f %*f %*d %*d %*d");
			start = false;
			Vec3f temp;
			temp.Set(x,y,z-10);
			line.push_back(temp);
		}
		fclose(fp);
	}

}


void Mesh::PaintVotes(int number)
{
	//printf("Begin PaintVotes\n");
	//	get_centerlines_voxels();
	float temp;
	Vec3f col;

	printf("Came here\n");
	for(unsigned int counter=0; counter<debug_vector.size(); counter++)
	{
		printf("Did I have anything in debug_vector[]?\n");
		glBegin(GL_LINES);
		glColor3f(1,0.5,0);
		glVertex3f(debug_vector[counter].start.x(),debug_vector[counter].start.y(),debug_vector[counter].start.z());
		glVertex3f(debug_vector[counter].destination.x(),debug_vector[counter].destination.y(),debug_vector[counter].destination.z());
		glEnd();
	}

	for(unsigned int counter=0; counter<debug_vector1.size(); counter++)
	{

		glBegin(GL_LINES);
		glColor3f(0,0,1);
		glVertex3f(debug_vector1[counter].start.x(),debug_vector1[counter].start.y(),debug_vector1[counter].start.z());
		glVertex3f(debug_vector1[counter].destination.x(),debug_vector1[counter].destination.y(),debug_vector1[counter].destination.z());
		glEnd();
	}

	for(unsigned int counter=0; counter<debug_vector2.size(); counter++)
	{

		glBegin(GL_LINES);
		glColor3f(0,1,0);
		glVertex3f(debug_vector2[counter].start.x(),debug_vector2[counter].start.y(),debug_vector2[counter].start.z());
		glVertex3f(debug_vector2[counter].destination.x(),debug_vector2[counter].destination.y(),debug_vector2[counter].destination.z());
		glEnd();
	}

	/*		
	glBegin(GL_POINTS);
	glColor3f(1,1,0);
	for(int counter=0; counter<debug_vector2.size(); counter++)
	{
	glVertex3f(debug_vector2[counter].x(),debug_vector2[counter].y(),debug_vector2[counter].z());
	}
	glEnd();
	*/
	glPointSize(2);
	glBegin(GL_POINTS);
	for(int cox =0; cox<maxx;cox++)
		for(int coy=0;coy<maxy;coy++)
			for(int coz=0;coz<maxz;coz++)
			{
#ifndef DISPFIT	
				temp = 	VOTES(cox,coy,coz);
				if(temp>=20+3*number)
				{
					//				printf("I came");
					col = getColorCode(temp,25+4*number,20+3*number,30+5*number);
					glColor3f(col.x(),col.y(),col.z());
					//glColor3f(0,1,0);
					glVertex3f(cox,coy,coz);
				}
#else
				temp = 	FIT(cox,coy,coz);
				if(temp<5 && VOTES(cox,coy,coz)>7)
				{
					//				printf("I came");
					col = getColorCode(temp,0.05,0,0.1);
					glColor3f(col.x(),col.y(),col.z());
					glVertex3f(cox,coy,coz);
				}
#endif 
			}
			glEnd();
			//printf("\n");

}

double cast_single_ray(Vec3f voxel, Vec3f direction, double max_distance,double max_radius)
{

	//	Iterator<Face*> *iter;
	//	Bag<Face*>*bf=m->getFaces();
	//	iter = bf->StartIteration();
	const double lambda = 1;
	Face * facep[2]={NULL,NULL};
	int maximum_steps = max_distance/lambda; // MAX STEPS

	//Array<Face *> arr(400);
	//	Array<Face *> arr1(400);
	//	Array<Face *> arr2(200);
	//	Array<Face *> arr3(200);
	//Array<Face *> temparray(50);
	//	set<Face*,lface> s1,s2,s3;
	//	vector <Face*> v1,v2,v3;
	//int number = 0;
	//double distance_opp =0;
	//int number123=3;

	//	scanf("%d",&number123);
	//number123 = 111111111;
	//number ++;
	/*(if(number % 1000==0)
	{
	printf("%d\n",number);
	//break;
	}*/
	//	printf("%d\n",number);
	//find the opposite
	Vec3f p = voxel;
	Vec3f n = Vec3f(0,0,0)-direction;//CHANGED!!!!
	n.Normalize();
	int steps =-1;
	Face *fminpos = NULL;
	double fmindist = 23432423;
	//	printf("Beginning search\n");
	Vec3f px = p - n;//*lambda
	while(steps < maximum_steps)
	{
		steps ++;
		px = px + n;//*lambda;// + lambda*steps*n;
		//px.Set((int)(px.x()+0.5),(int)(px.y()+0.5),(int)(px.z()+0.5));
		//	debug_vector2.push_back(px);
		const int w=1;
		for(int counterx = int(px.x()-w+0.5);counterx<int(px.x()+w+0.5);counterx++)
		{
			for(int countery = int(px.y()-w+0.5);countery<int(px.y()+w+0.5);countery++)
			{
				for(int counterz = int(px.z()-w+0.5);counterz<int(px.z()+w+0.5);counterz++)
				{
					Vec3f tempx = Vec3f(counterx,countery,counterz);
					//	Vec3f tempx = px;
					//	int count1=0;
					pair<map_type::const_iterator, map_type::const_iterator> pa = hash_for_centerline.equal_range(calcint(tempx));
					for (map_type::const_iterator i = pa.first; i != pa.second; ++i)
					{
						//		count1++;
						Face * fa = (*i).second;

						/*
						Vec3f nor = fa->getNormal();
						//	if(number == 3417 || number == 5323 || number == 9691 || number == 10869 || number == 14124)
						//		fa->setColor(Vec3f(1.0,0.0,1.0));


						nor.Normalize();
						if(n.Dot3(nor)>0)
						continue;
						*/
						//	if(fa->getColor()!=Vec3f(1,0,0)&&fa->getColor()!=Vec3f(0,0,0))
						//		fa->setColor(Vec3f(1,0,1));
						Vec3f a = getCenter(fa);
						double l = n.Dot3(a-p)/n.Dot3(n);

						if(l <0 )
						{
							continue;
							//printf("spurious result - check why this is happening\n");
						}

						double dist = (a-(p+l*n)).Length();
						if(dist <fmindist)
						{
							fmindist = dist;
							fminpos = fa;
						}

					}
					//	if(number==91)
					//	printf("count1 = %d\n",count1);
				}
			}
		}
		//	printf("Ended\n");
		if(fminpos!=NULL)
		{
			if(facep[0]==NULL)
			{
				/*twop twopoints;
				twopoints.start = voxel;
				twopoints.destination = fminpos->getCenter();
				debug_vector.push_back(twopoints);
				printf(" I did push into debug_vector\n");
				fminpos->setColor(Vec3f(1,0,0));
				*/

				facep[0]=fminpos;

				//	temparray.Clear();
				//	BFS(facep[0],5,&temparray);
				//fminpos = NULL;
				//fmindist = 1232131;
				break;
			}
			else
			{

				//int counter;
				//for(counter=0; counter<temparray.Count(); counter++)
				//{
				//	if(temparray[counter]==fminpos)
				//		break;
				//}
				//if(counter==temparray.Count())
				//{
				//	/*twop twopoints;
				//	twopoints.start = facep[0]->getCenter();
				//	twopoints.destination = fminpos->getCenter();
				//	debug_vector1.push_back(twopoints);
				//	*/
				//	facep[1]=fminpos;
				//	fminpos->setColor(Vec3f(0,0,1));
				//	break;
				//}
				fminpos= NULL;
				fmindist = 1223123;
			}
		}
		else
		{
			/*	twop twopoints;
			twopoints.start = voxel;
			twopoints.destination = px;
			debug_vector1.push_back(twopoints);*/
			//printf("Couldnt find\n");
		}

	}

	if(facep[0]!=NULL)
	{
		if (1)
			return fabs(facep[0]->getNormal().Dot3(n));
		else
			return 0;
		/*
		if(facep[0]->getNormal().Dot3(n)>0)
		return true;
		else
		return false;
		*/
	}
	else
	{
		/*twop twopoints;
		twopoints.start = voxel;
		twopoints.destination = p + lambda*maximum_steps*n;;
		debug_vector2.push_back(twopoints);
		*/
		return 0;
		//return false;
	}

}
double cast_rays(const Vec3f svoxel,const int n,const double max_distance,const double max_radius)
{
	//	printf("I did cast rays once here\n");
	const int theta_n = int(sqrt(float(n))+0.5); // round it
	const int phi_n = theta_n;//int(pow((float)n,(float)(1/2.0))+0.5);

	double theta_step=2*M_PI/theta_n;
	double phi_step = M_PI/phi_n;

	//int num_success = 0;
	double theta, phi;
	Vec3f direction;
	//int num_total = 0;
	double return_total=0;
	for(int theta_c=0; theta_c<theta_n; theta_c++)
	{
		theta = 0 + theta_c*theta_step;
		//#pragma omp parallel num_threads(2)
		for(int phi_c=0; phi_c<phi_n; phi_c++)
		{
			//	printf("Hi from %d\n",omp_get_thread_num());
			phi = 0 + phi_c*phi_step;
			direction.Set(cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi));
			return_total  += cast_single_ray(svoxel,direction,max_distance,max_radius);
			/*
			if(success)
			num_success++;
			*/
			//num_total++;
		}
	}

	//printf("num_total %d num_success %d\n",num_total,num_success);
	return return_total;//*1.0/num_total;
}
void scale_vector(double a, double b, double c)
{
	for(unsigned int counter=0; counter<debug_vector.size(); counter++)
	{
		twop temp = debug_vector[counter];
		temp.start.Set(temp.start[0]*a,temp.start[1]*b,temp.start[2]*c);
		debug_vector[counter]=temp;
	}
}
void get_centerlines_voxels()
{

	debug_vector.clear();
	debug_vector1.clear();
	debug_vector2.clear();
	srand(time(NULL));
	static double ratio_distance;

	//hash all the face centers rounded off to the nearest half-integers ( nice term :P )


	static bool once = false;
	static double dist = -1;
	static double distsum = 0;
	static int distnum = 0;
	double dist1,distmin;
	if(!once)
	{
		votes = new int[maxx*maxy*maxz];

		for(int counterx =0; counterx<maxx;counterx++)
		{
			for(int countery = 0; countery <maxy; countery++)
			{
				for(int counterz = 0; counterz < maxz; counterz++)
				{
					VOTES(counterx,countery,counterz)=0;	
				}
			}
		}
		Bag<Face*> * bf = m->getFaces();
		Iterator<Face *>*iter = bf->StartIteration();


		while(Face *f = iter->GetNext())
		{
			Vec3f center = getCenter(f);
			for(int counter =0; counter <3; counter++)
			{
				distsum += (center - (*f)[counter]->get()).Length();
				distnum++;
				if(dist < (center - (*f)[counter]->get()).Length())
				{
					dist = (center - (*f)[counter]->get()).Length();
				}
			}
		}
		distsum /= distnum;
		bf->EndIteration(iter);
		printf("maximum radius %lf\n",dist);
		//scale_vertices(0.3/dist,0.3/dist,0.3/dist);
		iter = bf->StartIteration();
		while(Face *f = iter->GetNext())
		{
			Vec3f center = getCenter(f);
			hash_for_centerline.insert(map_type::value_type(calcint(center),f));
		}
		bf->EndIteration(iter);
		//made the hash function! yo!


		iter = bf->StartIteration();
		dist1 = -1;distmin = 23453242;
		while(Face *f = iter->GetNext())
		{
			Vec3f center = getCenter(f);
			double max1 = -1;
			for(int counter =0; counter <3; counter++)
			{
				if(dist1 < (center - (*f)[counter]->get()).Length())
				{
					dist1 = (center - (*f)[counter]->get()).Length();
				}
				if(distmin > (center - (*f)[counter]->get()).Length())
				{
					distmin = (center - (*f)[counter]->get()).Length();
				}
				if(max1 < (center - (*f)[counter]->get()).Length())
				{
					max1 = (center - (*f)[counter]->get()).Length();
				}
			}
			if(max1 < 0.1)
			{
				//	f->setColor(Vec3f(1.0,0.0,1.0));
			}
		}
		printf("max radius %lf min radius %lf\n",dist1,distmin);
		ratio_distance = dist1/dist;
		//	scanf("%*d");
		bf->EndIteration(iter);
		once = true;
	}
	Vec3f rand_point;
	const int skipx=3;
	const int skipy=3;
	const int skipz=1;
	//int count_done = 0;
	printf("maxx maxy maxz %d %d %d\n",maxx,maxy,maxz);
	const int maxstepst = 35;

	for(int cox = MAX(minx-1,0); cox<maxx;cox+=skipx)
	{
		printf("cox %d  \n",cox);
		for(int coy = MAX(miny-1,0); coy<maxy;coy+=skipy)
		{
			printf("coy %d  \r",coy);
			fflush(stdout);
			//	
			//
//#pragma omp parallel for num_threads(2)
			for(int coz = MAX(minz-1,0); coz < maxz;coz+=skipz)
			{
				//rand_point.Set(rand()%(int(maxx*ratio_distance)),rand()%(int(maxy*ratio_distance)),rand()%(int(maxz*ratio_distance)));
				rand_point.Set(cox,coy,coz);
				int num = cast_rays(rand_point,100,maxstepst,2*distsum);
				VOTES(cox,coy,coz)=num;
				//return;
				//count_done++;
				//if(count_done>3)
				//	return;
			}
		}
	}
	//	scale_vector(1/ratio_distance,1/ratio_distance,1/ratio_distance);
	//	scale_vertices(1/ratio_distance,1/ratio_distance,1/ratio_distance);

}


void get_centerlines_cast_rays()
{
	Array<Face *> arr(50);
	Array<Face *> arr_new(50);	
	Bag<Face*>* bf = m->getFaces();
	Iterator<Face*>*iter = bf->StartIteration();
	int num_count =0;

	printf("maxx %d maxy %d maxz %d\n",maxx,maxy,maxz);

	votes = new int[maxx*maxy*maxz];

	fit_values = new float[maxx*maxy*maxz];
	for(int counter =0; counter<maxx*maxy*maxz;counter++)
	{
		fit_values[counter]=INT_MAX/2;
	}
	printf("Finished assigning fit values\n");
	float fit_value=100;
	while(Face *f = iter->GetNext())
	{
		num_count++;
		if(f->getArea()<0.05)
		{
			f->setColor(Vec3f(1,0,0));
			continue;
		}

		if(num_count %2000 ==0 )
			printf("%d\r",num_count);

		if(1)
		{
			arr.Clear();
			BFS_edge(f,2,&arr,true);
			//DEBUGP("Done BFS_edge\n");
			for(int counter=0; counter<arr.Count();counter++)
			{
				arr_new.Clear();
				BFS4(f,arr[counter],25,&arr_new,&arr);
				int size = arr_new.Count();
				for(int counter1=0; counter1<size; counter1++)
				{
					//	printf("Unmarking %d\n",counter1);
					arr_new[counter1]->unMark();
					arr_new[counter1]->setColor(Vec3f(0,1,0));
				}
				//DEBUGP("Done BFS4\n");
				Vec3f center = fit_circle(&arr_new,fit_value,20,2);
				//printf("fit = %f\n",fit_value);
				//center.Write();
				center = Vec3f(int(center.x()+0.5),int(center.y()+0.5),int(center.z()+0.5));
				Vec3f centerold = center;
				center.Set(MAX(MIN(center.x(),maxx-1),0),MAX(MIN(center.y(),maxy-1),0),MAX(MIN(center.z(),maxz-1),0));

				twop temptwop;
				if(center.Length()>5)
				{
					temptwop.start = arr_new[arr_new.Count()/2]->getCenter();
					temptwop.destination = center;
				}

				//	debug_vector.push_back(temptwop);
				int vw = 2;
				Vec3f minc = Vec3f(MAX(center[0]-vw,0),MAX(center[1]-vw,0),MAX(center[2]-vw,0));
				Vec3f maxc = Vec3f(MIN(center[0]+vw,maxx-1),MIN(center[1]+vw,maxy-1),MIN(center[2]+vw,maxz-1));
				for(int cox = minc[0];cox<maxc[0];cox++)
					for(int coy =minc[1];coy<maxc[1];coy++)
						for(int coz= minc[2];coz<maxc[2];coz++)
						{
							VOTES(cox,coy,coz)++;
							FIT(cox,coy,coz)=MIN(FIT(cox,coy,coz),fit_value);
						}
			}
		}
		if(num_count>5)
			break;
	}
	bf->EndIteration(iter);
}


void get_centerlines_normals()
{
	const int total_steps = 15;
	Array<Face *> arr(50);
	Array<Face *> arr_new(50);	
	Bag<Face*>* bf = m->getFaces();
	Iterator<Face*>*iter = bf->StartIteration();
	int num_count =0;

	printf("maxx %d maxy %d maxz %d\n",maxx,maxy,maxz);

	votes = new int[maxx*maxy*maxz];

	fit_values = new float[maxx*maxy*maxz];
	for(int counter =0; counter<maxx*maxy*maxz;counter++)
	{
		fit_values[counter]=INT_MAX/2;
	}
	printf("Finished assigning fit values\n");
	float fit_value=100;
	while(Face *f = iter->GetNext())
	{
		num_count++;
		if(num_count %2000 ==0 )
			printf("%d\r",num_count);

		if(1)
		{
			arr.Clear();
			BFS(f,1,&arr);
			//DEBUGP("Done BFS_edge\n");
			int size = arr.Count();
			Vec3f avgnormal = Vec3f(0,0,0);
			int avgn = 0;
			Vec3f tempnormal;
			for(int counter1=0; counter1<size; counter1++)
			{
				//	printf("Unmarking %d\n",counter1);
				arr[counter1]->unMark();
				tempnormal = arr[counter1]->getNormal();
				if(tempnormal.Length()<0.2)
					continue;
				tempnormal.Normalize();
				avgnormal = avgnormal+ tempnormal;
				avgn++;
			}
			avgnormal = avgnormal/avgn;

			avgnormal.Normalize();
			//Vec3f center = fit_circle(&arr_new,fit_value,20,2);
			Vec3f center = f->getCenter();

			for(int counter1=2; counter1 <= total_steps; counter1++)
			{
				//printf("fit = %f\n",fit_value);
				//center.Write();
				center=f->getCenter()-counter1*avgnormal;
				Vec3f centerold = center;

				center.Set(MAX(MIN(center.x(),maxx-1),0),MAX(MIN(center.y(),maxy-1),0),MAX(MIN(center.z(),maxz-1),0));

				if((center-centerold).Length()>1)
				{
					continue;
				}

				//twop temptwop;
				//if(center.Length()>5)
				//{
				//	temptwop.start = arr_new[arr_new.Count()/2]->getCenter();
				//	temptwop.destination = center;
				//}

				//debug_vector.push_back(temptwop);
				int vw = 1;
				Vec3f minc = Vec3f(MAX(center[0]-vw,0),MAX(center[1]-vw,0),MAX(center[2]-vw,0));
				Vec3f maxc = Vec3f(MIN(center[0]+vw,maxx-1),MIN(center[1]+vw,maxy-1),MIN(center[2]+vw,maxz-1));
				for(int cox = minc[0];cox<maxc[0];cox++)
					for(int coy =minc[1];coy<maxc[1];coy++)
						for(int coz= minc[2];coz<maxc[2];coz++)
						{
							VOTES(cox,coy,coz)++;
							FIT(cox,coy,coz)=MIN(FIT(cox,coy,coz),fit_value);
						}
			}

		}
		//if(num_count>5)
		//	break;
	}
	bf->EndIteration(iter);
}

void get_centerlines_sphere()
{
	Array<Face *> arr(50);
	Array<Face *> arr_new(50);	
	Bag<Face*>* bf = m->getFaces();
	Iterator<Face*>*iter = bf->StartIteration();
	int num_count =0;

	printf("maxx %d maxy %d maxz %d\n",maxx,maxy,maxz);

	votes = new int[maxx*maxy*maxz];

	fit_values = new float[maxx*maxy*maxz];
	for(int counter =0; counter<maxx*maxy*maxz;counter++)
	{
		fit_values[counter]=INT_MAX/2;
	}
	printf("Finished assigning fit values\n");
	float fit_value=100;
	while(Face *f = iter->GetNext())
	{
		num_count++;
		if(num_count %2000 ==0 )
			printf("%d\r",num_count);

		if(1)
		{
			arr.Clear();
			BFS(f,8,&arr);
			//DEBUGP("Done BFS\n");
			for(int counter=0; counter<arr.Count();counter++)
				arr[counter]->setColor(Vec3f(1,0,0));
			arr_new.Clear();


			Vec3f center = fit_circle(&arr,fit_value,100,2);
			//printf("fit = %f\n",fit_value);
			//center.Write();
			center = Vec3f(int(center.x()+0.5),int(center.y()+0.5),int(center.z()+0.5));
			Vec3f centerold = center;
			center.Set(MAX(MIN(center.x(),maxx-1),0),MAX(MIN(center.y(),maxy-1),0),MAX(MIN(center.z(),maxz-1),0));

			twop temptwop;
			if(center.Length()>5)
			{
				temptwop.start = f->getCenter();
				temptwop.destination = center;
			}

			//	debug_vector.push_back(temptwop);
			int vw = 2;
			Vec3f minc = Vec3f(MAX(center[0]-vw,0),MAX(center[1]-vw,0),MAX(center[2]-vw,0));
			Vec3f maxc = Vec3f(MIN(center[0]+vw,maxx-1),MIN(center[1]+vw,maxy-1),MIN(center[2]+vw,maxz-1));
			for(int cox = minc[0];cox<maxc[0];cox++)
				for(int coy =minc[1];coy<maxc[1];coy++)
					for(int coz= minc[2];coz<maxc[2];coz++)
					{
						VOTES(cox,coy,coz)++;
						FIT(cox,coy,coz)=MIN(FIT(cox,coy,coz),fit_value);
					}

		}
		//	if(num_count>10)
		//		break;
	}
	bf->EndIteration(iter);
}

void write_votes_to_file( char *filename)
{
	FILE *fp = fopen(filename,"w");
	for(int cox =0; cox <maxx; cox++)
	{
		for(int coy = 0; coy<maxy; coy++)
		{
			for(int coz =0; coz <maxz; coz++)
			{
				if(VOTES(cox,coy,coz)>0)
				{
					int cx = 0;
					int cy = 0;
					for(cx = -1; cx<=1 ; cx++)
					{
						for(cy = -1; cy<=1; cy++)
						{
							if(cox+cx < maxx && cox+cx >=0 && coy+cy<maxy && coy+cy>=0)
								fprintf(fp,"%d %d %d %d\n",cox+cx,coy+cy,coz,VOTES(cox,coy,coz));
						}
					}
				}
			}
		}
	}
	fclose(fp);
}

void load_votes_from_file(char *filename)
{
	FILE *fp = fopen(filename,"r");
	int a[4];
	//float b;

	if(maxx == 0 || maxy ==0 || maxz==0)
	{
		printf("Load maxx maxy and maxz first and then call load_votes_from_file \n");
	}
	votes = new int[maxx*maxy*maxz];

	if (votes == NULL)
	{
		printf("Couldnt allocate memory for votes\n");
	}
	
	
	fit_values = new float[maxx*maxy*maxz];

	for(int cox =0; cox <maxx; cox++)
	{
		for(int coy = 0; coy<maxy; coy++)
		{
			for(int coz =0; coz <maxz; coz++)
			{
				if( fscanf(fp,"%d %d %d %d",&a[0],&a[1],&a[2],&a[3]) == EOF )
          {
          cerr << "EOF encountered in fscanf" << endl;
          }
				if(feof(fp))
					goto label;
			//	printf("%d %d %d %d\n",a[0],a[1],a[2],a[3]);
				VOTES(a[0],a[1],a[2])=a[3];
			//	FIT(a[0],a[1],a[2])=b;
			}
		}
	}
label:
	fclose(fp);
}

int main(int argc, char* argv[])
{

	current_scaling = Vec3f(1.0,1.0,1.0);
	int dec ;
	dec =100000;

	//    load_points("
	m = new Mesh();
	//	int c;

	for(int x = 0;x <1030;x++)
		for(int y = 0;y <1030;y++)
			for(int z = 0;z <82;z++)
				pmatrix[x][y][z]=0; 


	//char filename[1024] ;
	// "E:\\Arun\\My matlab codes\\deconvolved-lower3-higher8-window2";
	//char filename_data[1024];
	//char picfilename [1024];
	//strcpy(filename_data,"E:\\Arun\\My matlab codes\\phantom-binary\\gaussian_noise\\stdeviation_0");
	//strcat(filename_data,".npts");
	//strcpy(picfilename,"E:\\Arun\\My matlab codes\\phantom-binary\\gaussian_noise\\stdeviation_0");
	//strcat(picfilename,".pic");
	//m->set_pic_name(picfilename);

	//load_points_without_normal("deconvolved-good.npts");
	if(argc < 3)
	{
		printf("Usage centerline inputfile outputfile\n");
		return 1;
	}
	
	load_points_without_normal(argv[1]);

//	m->set_pic_name("D:\\MSA paper\\vessel.pic");
	//	read_from_file("../ascii.txt");
	printf("done reading from file\n");
	//	make_hypothesis_model();




	//	rank70Filter();


	marchtetra();

	//  LoadOBJ(argv[1]);
	//	LoadOBJ("thinned-without-color-for-confidence-measurement.obj");
	//	LoadOBJ("decimated-full-112.obj");
	//	LoadOBJ("try.obj");
	//	LoadOBJ("try-fullvessel.obj");
	//	LoadOBJ("deconvolved-lower3-higher10-window2.obj");
	//LoadOBJ("centerline-testing.obj");
	//	get_distance_statistics();
	//	check_consistent();
	//	delete_small_components();
	//	printf("Going for the second time to delete_small_components\n");
	//	delete_small_components();

	check_consistent();
	//scanf("%*d");

	scale_vertices(current_scaling.x(),current_scaling.y(),2*current_scaling.z());
	maxz = maxz*2;

	for(int counter =0; counter <5; counter++)
	{
		//try both laplacian and volume preserved smoothing
		//laplacian smoothing is good enough for centerline extraction since we want to obtain only the topology
		// it saves time
		if(0)
			SmoothSurface_complex(1);
		else
			SmoothSurface();
		printf("Initial smoothing %d time(s)\n", counter+1);
		if(0)
		{
			char buff[] = "latest-thinning-results-smoothed";
			char buffnew[1024];
			sprintf(buffnew,"%s-%d.obj",buff,counter/3);
			writeOBJ(buffnew);
		}
	}

	//find_min_median_max();	
	dec = 8000;
	printf("Num vertices (before decimating) %d\n",m->numVertices());
	printf("Num Faces (before decimating) %d\n",m->numFaces());
	printf("dec?");
	dec = 0;

	for(int counter=0; counter<1; counter++)
	{
		dec = 3*m->numFaces()/8;
		//scanf("%d",&dec);
		//	dec = 1000000;
		Decimate(dec);
	}
	check_consistent();


	//CENTERLINE EXTRACTION -uncomment load_votes_from_file() to load from a saved file.
	get_centerlines_voxels();
	//load_votes_from_file("half_vessel_votes5.txt");
	write_votes_to_file(argv[2]);

	char tif_filename[1024];
	strcpy(tif_filename,argv[2]);
	strcpy(&tif_filename[strlen(argv[2])-3],"tif");
	printf("tif_filename %s inputfilename %s outputfilename %s\n",tif_filename,argv[1],argv[2]);

//	char votesbuff[2000];
//	sprintf(votesbuff,"./votestotiff %s %s",argv[2],tif_filename);
//	printf("to run |%s|\n",votesbuff);
//	system(votesbuff);
	

	return 0;

	printf("Num vertices %d\n",m->numVertices());
	printf("Num Faces %d\n",m->numFaces());
	printf("returned_count %d, non_returned_count %d\n",returned_count,non_returned_count);
	// facebag->EndIteration(iter);
	ArgParser *ap = new ArgParser();
	//scanf("%*d");
	printf("Finished Initializing ArgParser .....\n");
	//scanf("%*d");
	Radiosity* r= new Radiosity(m,ap);
	printf("Finished Initializing Radiosity .....\n");
	//scanf("%*d");
	//    glBegin(GL_LINES);
	GLCanvas *g;
	//BoundingBox *b = m->getBoundingBox();
	//      m->PaintWireframe();
	g->initialize(ap,r);

	//	printf("hi\n");

	//m->PaintWireframe();
	// b->paint();
	return 0;
}
