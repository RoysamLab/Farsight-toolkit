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

#include <stdio.h>
#include <assert.h>
#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif 
#include "mesh.h"
#include "edge.h"
#include "vertex.h"
#include "face.h"
#include "glCanvas.h"
#include "vertex_parent.h"

#define INITIAL_VERTEX 10000
#define INITIAL_EDGE 10000
#define INITIAL_FACE 10000

// =======================================================================
// CONSTRUCTORS & DESTRUCTORS
// =======================================================================

Mesh::Mesh() {
	
  vertices = new Array<Vertex*>(INITIAL_VERTEX);
  edges = new Bag<Edge*>(INITIAL_EDGE,Edge::extract_func);
  faces = new Bag<Face*>(INITIAL_FACE,Face::extract_func);
  vertex_parents = new Bag<VertexParent*>(INITIAL_VERTEX, VertexParent::extract_func);
  bbox = NULL;
}

Mesh::~Mesh() {
  delete vertices;  
  vertices = NULL;
  delete edges;
  edges = NULL;
  delete faces;
  faces = NULL;
  delete bbox;
  bbox = NULL;
 glDeleteTextures(rdepth,texturexy);  
}



void Mesh::set_pic_name(char *a)
{
	strcpy(pic_name,a);
}
void Mesh::texture_list_init()
{
	HandleGLError();
	printf("I am in texture_list_init\n");
	//scanf("%*d");
	FILE * fpi = fopen(pic_name,"rb");
	unsigned char temp[77];
	if( fread((void*)temp,sizeof(unsigned char),76,fpi) != 76 )
    {
    cerr << "Less than 76 elements read by fread" << endl;
    }
	//printf("%d %d %d\n",CFH(temp[0],0)+CFH(temp[1],2),CFH(temp[2],0)+CFH(temp[3],2),CFH(temp[4],0)+CFH(temp[5],2));
	rwidth = CFH(temp[0],0)+CFH(temp[1],2);
	rlength = CFH(temp[2],0)+CFH(temp[3],2);
	rdepth = CFH(temp[4],0)+CFH(temp[5],2);
	int twopow = 0;
	int rdepth1 = rdepth-1;
	while(rdepth1)
	{
		rdepth1= rdepth1>>1;
		twopow++;
	}
	//twopow--;
	rdepth1 = 1<<twopow;
	int npixels = rwidth*rlength;
	//int mulwd = rwidth*rdepth1;
	//int mulld = rlength*rdepth1;
	//int w = rwidth;
	
	unsigned char *rasterxy;// = (unsigned char*)malloc(npixels*(rdepth1)*sizeof(unsigned char));
	//unsigned char **rasteryz;// = (unsigned char*)malloc(npixels*(rdepth1)*sizeof(unsigned char));	
	//unsigned char **rasterxz;// = (unsigned char*)malloc(npixels*(rdepth1)*sizeof(unsigned char));
	rasterxy = (unsigned char*) malloc(rdepth*npixels*sizeof(unsigned char));
	/*
	for(int counter=0; counter<rdepth; counter++)
	{
		rasterxy[counter]= (unsigned char*)malloc(npixels*sizeof(unsigned char));
	}
	printf("allocated memory for rasterxy\n");
	
	rasteryz = (unsigned char**) malloc(rwidth*sizeof(unsigned char *));
	for(int counter=0; counter<rwidth; counter++)
	{
		rasteryz[counter]= (unsigned char*)malloc(mulld*sizeof(unsigned char));
	}
	printf("allocated memory for rasteryz\n");
	rasterxz = (unsigned char**) malloc(rlength*sizeof(unsigned char *));
	for(int counter=0; counter<rlength; counter++)
	{
		rasterxz[counter]= (unsigned char*)malloc(mulwd*sizeof(unsigned char));
	}
	printf("allocated memory for rasterxz\n");
	*/
	//allocate_memory(rasterxy,rdepth,rlength*rwidth);
	//allocate_memory(rasteryz,rwidth,rlength*rdepth);
	//allocate_memory(rasterxz,rlength,rwidth*rdepth);

	if(rasterxy==NULL)
	{
		printf("memory problem: couldn't allocate enough memory\n");
		_exit(0);
	}
	//printf("%d %d %d\n",rasterxy, rasteryz, rasterxz);
	printf("I'm going to read the memory block now rdepth = %d rheight = %d rwidth = %d\n", rdepth, rlength, rwidth);
//	for(int counter=0; counter<rdepth;counter++)
//		fread(rasterxy[counter],sizeof(unsigned char),npixels,fpi);
	printf("I read %d\n",(int)fread(rasterxy,sizeof(unsigned char),npixels*rdepth,fpi));
	//rdepth+=1;//CHANGED
	printf("done reading %d %d %d %d %zu\n",rwidth, rlength,rdepth,(1<<twopow),sizeof(unsigned char));
	printf("I'm encoding the data now\n");
	//double max1=-1,min1=23423423,mean1=0;
	//int num1;
	//int pc = 0;
	//unsigned char test;
	//rdepth = rdepth1;
	printf("%d rdepth\n",rdepth);
	/*
	for(int cx = 0; cx<rwidth; cx++)
	{
		printf("%d \r",cx);
		for(int cy =0; cy<rlength; cy++)
			for(int cz = 0; cz<rdepth;cz++) 
			{
				//test = rasterxy[cz*npixels+cy*rwidth+cx];
				rasteryz[cx][cy*rdepth+cz]=rasterxy[cz][cy*rwidth+cx];
				rasterxz[cy][cx*rdepth+cz]=rasterxy[cz][cy*rwidth+cx];
				//rasteryz[(pc - pc%(mulld))+(pc%(mulld))/rlength+pc%rdepth]=test;
				//rasterxz[(pc-pc%(mulwd))+(pc%(mulwd))/rwidth+pc%rdepth]=test;
				//if(max1<test)
				//	max1=test;
				//if(min1>test)
				//	min1=test;
				//mean1=mean1+test;
				//raster[cz*npixels+cy*rwidth+cx]=255-raster[cz*npixels+cy*rwidth+cx];
				pc ++;
			}
	}*/
	//printf("min mean max %lf %lf %lf\n",min1,mean1/num1,max1);
	//texturexy=NULL;
	//textureyz=NULL;
	//texturexz=NULL;
	texturexy = new GLuint[rdepth];
	glDeleteTextures(rdepth,texturexy);
	glGenTextures(rdepth,texturexy);
	HandleGLError();
	printf("Finished generating textures.. now loading data\n");
	//textureyz = texturexy +rdepth;
	//texturexz = texturexy +10+rwidth;
	/*glGenTextures(1,textureyz);
	glGenTextures(1,texturexz);*/
	//printf("pointers %d %d %d\n",texturexy,textureyz,texturexz);
	list = glGenLists(1);
	glEnable(GL_TEXTURE_2D);
	//glEnable(GL_TEXTURE_3D);
	HandleGLError();
	glPixelTransferf(GL_RED_SCALE,5.0);
	for(int counter=0; counter<rdepth;counter++)
	{
	//	printf("%u ",texturexy[counter]);
//	glBindTexture(GL_TEXTURE_3D,texture[0]);
	glBindTexture(GL_TEXTURE_2D,texturexy[counter]);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	//PFNGLTEXIMAGE3DPROC glTexImage3D;
	//glTexImage3D = (PFNGLTEXIMAGE3DPROC) wglGetProcAddress("glTexImage3D");
	//printf("hi");
	printf(".");
	//HandleGLError();
	glTexImage2D(GL_TEXTURE_2D,0,GL_INTENSITY,rwidth,rlength,0,GL_RED,GL_UNSIGNED_BYTE,rasterxy+npixels*counter);
	//glTexImage3D(GL_TEXTURE_3D, 0, 1,rwidth,rlength,1<<twopow,0,GL_RED,GL_UNSIGNED_BYTE,raster);
	}
	printf("\n");
	/*
	for(int counter=0; counter<rwidth;counter++)
	{
		printf("si\n");
		printf("%u \n",textureyz[counter]);
		glBindTexture(GL_TEXTURE_2D,textureyz[counter]);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
		printf("si start\n");
		glTexImage2D(GL_TEXTURE_2D,0,GL_INTENSITY,rdepth,rlength,0,GL_RED,GL_UNSIGNED_BYTE,rasteryz[counter]);
		printf("si1\n");
		HandleGLError();
	}
	for(int counter=0; counter<rlength;counter++)
	{
		printf("qi\n");
		glBindTexture(GL_TEXTURE_2D,texturexz[counter]);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D,0,GL_INTENSITY,rdepth,rwidth,0,GL_RED,GL_UNSIGNED_BYTE,rasterxz[counter]);
	}
	printf("hi1\n");
	*/
	HandleGLError();
		GLfloat f;
	glGetFloatv(GL_RED_SCALE,&f);
	printf("GL_RED SCALE %lf\n",f);
	rdepth1= 1<<twopow;
     /*glNewList(list,GL_COMPILE_AND_EXECUTE);
     glBindTexture(GL_TEXTURE_3D, texture);
     glBegin(GL_QUADS);
     glTexCoord3d(0,0,0);glVertex3f(0,0,0);
     glTexCoord3d(1,0,0);glVertex3f(rwidth,0,0);
     glTexCoord3d(0,1,0);glVertex3f(0,rlength,0);
     glTexCoord3d(0,0,1);glVertex3f(0,0,rdepth1);
     glTexCoord3d(1,1,0);glVertex3f(rwidth,rlength,0);
     glTexCoord3d(1,0,1);glVertex3f(rwidth,0,rdepth1);
     glTexCoord3d(0,1,1);glVertex3f(0,rlength,rdepth1);
     glTexCoord3d(1,1,1);glVertex3f(rwidth,rlength,rdepth1);
     glEnd();
     glEndList();*/
     printf("Exiting texture_list_init \n");	
	 free(rasterxy);
	 //free(rasterxz);
	 //free(rasteryz);
	//double si = 10;

}
// =======================================================================
// MODIFIERS:   ADD & REMOVE
// =======================================================================

Vertex* Mesh::addVertex(const Vec3f &position) {
  int index = numVertices();
  //int tempx = int(10.0*(position.x()-int(position.x()))+0.5)%10;
//  int tempy = int(10.0*(position.y()-int(position.y()))+0.5)%10;
//  int tempz = int(10.0*(position.z()-int(position.z()))+0.5)%10;
//  
  //if(!(tempx == 1 || tempx==9 || tempx==0) || !(tempy == 1 || tempy==9 || tempy==0) || !(tempz == 1 || tempz==9 ||tempz==0))
//  {
//		printf("I got %lf %lf %lf %d %d %d\n",position.x(),position.y(),position.z(),tempx,tempy,tempz);
//		
//  }
  Vertex *v = new Vertex(index, position);
  vertices->Add(v);
  if (bbox == NULL) 
    bbox = new BoundingBox(position,position);
  else 
    bbox->Extend(position);
  return v;
}

void Mesh::addFace(Vertex *a, Vertex *b, Vertex *c, const Vec3f &color, const Vec3f &emit) {

	
  // create the face
  Face *f = new Face(color,emit);
assert(f);
  // create the edges
  Edge *ea = new Edge(a,f);
  Edge *eb = new Edge(b,f);
  Edge *ec = new Edge(c,f);
  //Edge *ed = new Edge(d,f);
  assert(ea);
  assert(eb);
  assert(ec);
  // point the face to one of its edges
  f->setEdge(ea);

  // connect the edges to each other
  ea->setNext(eb);
  eb->setNext(ec);
  ec->setNext(ea);
//  ed->setNext(ea);

  // add them to the master list
  edges->Add(ea);
  edges->Add(eb);
  edges->Add(ec);
  //edges->Add(ed);

  // connect up with opposite edges (if they exist)
  Edge *ea_op = getEdge((*ea)[1],(*ea)[0]);
  Edge *eb_op = getEdge((*eb)[1],(*eb)[0]);
  Edge *ec_op = getEdge((*ec)[1],(*ec)[0]);  
  //Edge *ed_op = getEdge((*ed)[1],(*ed)[0]);  

  if (ea_op != NULL) { ea_op->setOpposite(ea); }
  if (eb_op != NULL) { eb_op->setOpposite(eb); }
  if (ec_op != NULL) { ec_op->setOpposite(ec); }
//  if (ed_op != NULL) { ed_op->setOpposite(ed); }

  // add the face to the master list
  faces->Add(f); 
}
Face * Mesh::addFace2(Vertex *a, Vertex *b, Vertex *c, const Vec3f &color, const Vec3f &emit) {

  // create the face
  Face *f = new Face(color,emit);

  // create the edges
  Edge *ea = new Edge(a,f);
  Edge *eb = new Edge(b,f);
  Edge *ec = new Edge(c,f);

  //Edge *ed = new Edge(d,f);

  // point the face to one of its edges
  f->setEdge(ea);

  // connect the edges to each other
  ea->setNext(eb);
  eb->setNext(ec);
  ec->setNext(ea);
//  ed->setNext(ea);

  // add them to the master list
  edges->Add(ea);
  edges->Add(eb);
  edges->Add(ec);
  //edges->Add(ed);

  // connect up with opposite edges (if they exist)
 
  Edge *ea_op = getEdge((*ea)[1],(*ea)[0]);
  Edge *eb_op = getEdge((*eb)[1],(*eb)[0]);
  Edge *ec_op = getEdge((*ec)[1],(*ec)[0]);  
  
  //Edge *ed_op = getEdge((*ed)[1],(*ed)[0]);  

	
  if (ea_op != NULL) { ea_op->setOpposite(ea); }
  if (eb_op != NULL) { eb_op->setOpposite(eb); }
 if (ec_op != NULL) { ec_op->setOpposite(ec); }
//  for(;;);
//  if (ed_op != NULL) { ed_op->setOpposite(ed); }

  // add the face to the master list
  faces->Add(f); 
  return f;
}
void Mesh::removeFace(Face *f) {

//printf("Entered here\n");
  Edge *ea = f->getEdge();
  Edge *eb = ea->getNext();
  Edge *ec = eb->getNext();
  
  //Edge *ed = ec->getNext();


  assert (ec->getNext() == ea);

  // remove elements from master lists
  
  edges->Remove(ea);
  edges->Remove(eb);
  edges->Remove(ec);
  
//  edges->Remove(ed);
 //printf("just before faces->Remove(f)\n");

	//scanf("%*d");
	faces->Remove(f);
	//scanf("%*d");
	
  // clean up memory
  
  //ea->clearOpposite();
  //eb->clearOpposite();
 // ec->clearOpposite();
  
  delete ea;
  delete eb;
  delete ec;
  
 // delete ed;
  delete f;
 // for(int i =0;i<50000;i++)
  //		  printf("testing %d\r",i);
  //printf("out of remove face\n");
}

void Mesh::removeFace2(Face *f) {

  Edge *ea = f->getEdge();
  Edge *eb = ea->getNext();
  Edge *ec = eb->getNext();
  
  //Edge *ed = ec->getNext();

  assert (ec->getNext() == ea);

  // remove elements from master lists
//  printf("%d %d %d\n",ea,eb,ec);
  int unused = scanf("%*d");
  unused++;
  edges->Remove(ea);
  unused = scanf("%*d");
  unused++;
  edges->Remove(eb);
  unused = scanf("%*d");
  unused++;
  edges->Remove(ec);
  unused = scanf("%*d");
  unused++;
//  edges->Remove(ed);
 //printf("just before faces->Remove(f)\n");

	
	faces->Remove(f);
	
	
  // clean up memory
  delete ea;
  delete eb;
  delete ec;
 // delete ed;
  delete f;
 // for(int i =0;i<50000;i++)
  //		  printf("testing %d\r",i);
  //printf("out of remove face\n");
}
Edge* Mesh::getEdge(Vertex *a, Vertex *b) const {
  assert (edges != NULL);
  Edge *answer = edges->Get(a->getIndex(),b->getIndex());
  return answer;
}

Vertex* Mesh::getChildVertex(Vertex *p1, Vertex *p2) const {
  VertexParent *vp = vertex_parents->GetReorder(p1->getIndex(), p2->getIndex());
  if (vp == NULL) return NULL;
  return vp->get();
}

void Mesh::setParentsChild(Vertex *p1, Vertex *p2, Vertex *child) {
  vertex_parents->Add(new VertexParent(p1,p2,child));
}

//
// =======================================================================
// the load function parses very simple .obj files
// =======================================================================

void Mesh::Load(const char*input_file){}
//void Mesh::Load(const char *input_file) {
//  
//  FILE *objfile = fopen(input_file,"r");
//  if (objfile == NULL) {
//    printf ("ERROR! CANNOT OPEN '%s'\n",input_file);
//    return;
//  }
//
//  char line[200];
//  char token[100];
//  char atoken[100];
//  char btoken[100];
//  char ctoken[100];
//  char dtoken[100];
//  char etoken[100];
//  float x,y,z;
//  int a,b,c,d,e;
//  
//  int index = 0;
//  int vert_count = 0;
//  int vert_index = 1;
//
//  Vec3f color = Vec3f(1,1,1);
//  Vec3f emit = Vec3f(0,0,0);
//  
//  while (fgets(line, 200, objfile)) {   
//    
//    if (line[strlen(line)-2] == '\\') {
//      fgets(token, 100, objfile);	
//      int tmp = strlen(line)-2;
//      strncpy(&line[tmp],token,100);
//    }
//    int token_count = sscanf (line, "%s\n",token);
//    if (token_count == -1) continue;
//    a = b = c = d = e = -1;
//    if (!strcmp(token,"usemtl") ||
//	!strcmp(token,"g")) {
//      vert_index = 1; //vert_count + 1;
//      index++;
//    } else if (!strcmp(token,"v")) {
//      vert_count++;
//      sscanf (line, "%s %f %f %f\n",token,&x,&y,&z);
//      addVertex(Vec3f(x,y,z));
//    } else if (!strcmp(token,"f")) {
//      int num = sscanf (line, "%s %s %s %s %s %s\n",token,
//			atoken,btoken,ctoken,dtoken,etoken);
//      assert (num == 5);
//      sscanf (atoken,"%d",&a);
//      sscanf (btoken,"%d",&b);
//      sscanf (ctoken,"%d",&c);
//      sscanf (dtoken,"%d",&d);
//      a -= vert_index;
//      b -= vert_index;
//      c -= vert_index;
//      d -= vert_index;
//      assert (a >= 0 && a < numVertices());
//      assert (b >= 0 && b < numVertices());
//      assert (c >= 0 && c < numVertices());
//      assert (d >= 0 && d < numVertices());
//      addFace(getVertex(a),getVertex(b),getVertex(c),getVertex(d),color,emit);
//    } else if (!strcmp(token,"m")) {
//      // this is not standard .obj format!!
//      // materials
//      float r,g,b;
//      int num = sscanf (line, "%s %f %f %f", token,&r,&g,&b);
//      assert (num == 4);
//      color = Vec3f(r,g,b);
//    } else if (!strcmp(token,"l")) {
//      // this is not standard .obj format!!
//      // lights
//      float r,g,b;
//      int num = sscanf (line, "%s %f %f %f", token,&r,&g,&b);
//      assert (num == 4);
//      emit = Vec3f(r,g,b);
//    } else if (!strcmp(token,"vt")) {
//    } else if (!strcmp(token,"vn")) {
//    } else if (token[0] == '#') {
//    } else {
//      printf ("LINE: '%s'",line);
//    }
//  }
//}

// =======================================================================
// PAINT
// =======================================================================

void Mesh::PaintWireframe() {

  glDisable(GL_LIGHTING);

  // draw all the interior edges
  glLineWidth(1);
  glColor3f(1,1,1);
  glBegin (GL_LINES);
  Iterator<Edge*> *iter = edges->StartIteration();
  while (Edge *e = iter->GetNext()) {
    if (e->getOpposite() == NULL) continue;
    Vec3f a = (*e)[0]->get();
    Vec3f b = (*e)[1]->get();
    glVertex3f(a.x(),a.y(),a.z());
    glVertex3f(b.x(),b.y(),b.z());
  }
  edges->EndIteration(iter);
  glEnd();
  
  // draw all the boundary edges
//  glLineWidth(3);
//  glColor3f(1,0,0);
//  glBegin (GL_LINES);
//  iter = edges->StartIteration();
//  while (Edge *e = iter->GetNext()) {
//    if (e->getOpposite() != NULL) continue;
//    Vec3f a = (*e)[0]->get();
//    Vec3f b = (*e)[1]->get();
//    glVertex3f(a.x(),a.y(),a.z());
//    glVertex3f(b.x(),b.y(),b.z());
//  }
//  edges->EndIteration(iter);
//  glEnd();
//  
  glEnable(GL_LIGHTING);
  
  HandleGLError(); 
}

void Mesh::PaintPoints() {

  glDisable(GL_LIGHTING);
  float a,b;
  glGetFloatv(GL_POINT_SIZE_RANGE,&a);
  glGetFloatv(GL_POINT_SIZE_GRANULARITY,&b);
  printf("Point size range allowed %lf , min and max being %lf %lf\n", a-b,b,a);

  glLineWidth(1);
  glColor3f(1,0,0);
  glBegin (GL_POINTS);
	int counter =0;
  while (counter < vertices->Count()) {
		Vertex *v = (*vertices)[counter];
		counter ++;
    Vec3f a = v->get();
    glVertex3f(a.x(),a.y(),a.z());
  }

  glEnd();
  
  glEnable(GL_LIGHTING);
  
  HandleGLError(); 
}

// =================================================================
// SUBDIVISION
// =================================================================

Vertex* Mesh::AddEdgeVertex(Vertex *a, Vertex *b) {
  Vertex *v = getChildVertex(a,b);
  if (v != NULL) return v;
  Vec3f pos = 0.5*a->get() + 0.5*b->get();
  v = addVertex(pos);
  setParentsChild(a,b,v);
  return v;
}

Vertex* Mesh::AddMidVertex(Vertex *a, Vertex *b, Vertex *c) {
  Vec3f pos = 0.3333333*a->get() + 0.3333333*b->get() + 0.333333*c->get();
  Vertex *v = addVertex(pos);
  return v;
}

void Mesh::Subdivision() {
  printf ("Subdivide the mesh!\n");
  int i;

  Array<Face*> todo = Array<Face*>(numFaces());

  Iterator<Face*> *iter = faces->StartIteration();
  while (Face *f = iter->GetNext()) {
    todo.Add(f);
  }
  faces->EndIteration(iter);
 printf("%d faces to be split\n",todo.Count());
  for (i = 0; i <todo.Count(); i++) {
    Face *f = todo[i];

    Vertex *a = (*f)[0];
    Vertex *b = (*f)[1];
    Vertex *c = (*f)[2];
   // Vertex *d = (*f)[3];

    // add new vertices on the edges
    Vertex *ab = AddEdgeVertex(a,b);
    Vertex *bc = AddEdgeVertex(b,c);
    Vertex *ca = AddEdgeVertex(c,a);
   // Vertex *da = AddEdgeVertex(d,a);

    // add new point in the middle of the patch
    Vertex *mid = AddMidVertex(a,b,c);

    assert (getEdge(a,b) != NULL);
    assert (getEdge(b,c) != NULL);
    assert (getEdge(c,a) != NULL);
//    assert (getEdge(d,a) != NULL);

    // copy the color and emission from the old patch to the new
    Vec3f color = f->getColor();
    Vec3f emit = f->getEmit();
    removeFace(f);

    // create the new faces
    addFace(a,ab,mid,color,emit);
    addFace(b,mid,ab,color,emit);
    addFace(b,bc,mid,color,emit);
    addFace(mid,bc,c,color,emit);
    addFace(c,ca,mid,color,emit);
    addFace(mid,ca,a,color,emit);
   // addFace(d,da,mid,cd,color,emit);

//    assert (getEdge(a,ab) != NULL);
//    assert (getEdge(ab,b) != NULL);
//    assert (getEdge(b,bc) != NULL);
//    assert (getEdge(bc,c) != NULL);
//    assert (getEdge(c,ca) != NULL);
//    assert (getEdge(ca,a) != NULL);

   // assert (getEdge(d,da) != NULL);
   // assert (getEdge(da,a) != NULL);
  }
}

// =================================================================
