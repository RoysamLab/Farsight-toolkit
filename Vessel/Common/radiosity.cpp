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

#include "radiosity.h"
#include "mesh.h"
#include "face.h"
#include "glCanvas.h"
#include "draw_sphere.h"
#include "array.h"
#include <queue>
#include <algorithm>
#include <string.h>
#include "matrix.h"
using namespace std;
#define BUFSIZE 1000000

#ifndef M_PI
#define M_PI (2*acos(double(0)))
#endif

int GLCanvas::curr_type; 
GLfloat GLCanvas::xtrans,GLCanvas::ytrans,GLCanvas::ztrans;
Face *faces_picked[3];
int num_picked[3];
char curr_annotation[512];
int picked_faces  =0;

#define SKIP 2
// ================================================================
// CONSTRUCTOR & DESTRUCTOR
// ================================================================
void BFS_debug (Face *f, int depth, Array<Face*> *array)
{
	
	//printf("I am in BFS_debug and I got Face %p depth %d array_faces %p\n",f,depth,array);
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
		//if(q.size()<10)
		//printf("queue %d\t",q.size());
		temp = q.front();
		int index = qindex.front();
		q.pop();
		qindex.pop();
		if(index==0)
		{
			continue;
		}
		Edge *e = temp->getEdge();
		if(e==NULL)
			printf("e is null.. :((\n");
		//Edge * temp1 = NULL;
		for (int counter =0; counter < 3; counter ++)
		{
			Edge *e1 = e->getOpposite();
			if(e1!= NULL)
			{
				t = e1->getFace();
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
			if(e==NULL)
				printf("e is null now!\n");
		}
	}
//	printf("I'm at the end of debug_BFS\n");
}
void draw_transparent_sphere(GLfloat);
GLuint selectBuf[BUFSIZE];
Vec3f LAST_MOUSE;
Vec3f CENTER;
Face *PICKED = NULL;

Radiosity::Radiosity(Mesh *m, ArgParser *a) {
  mesh = m;
  args = a;
  Initialize();
  volume_rendering_setup_complete=false;
}

Radiosity::~Radiosity() {
	
  Cleanup();

}

void Radiosity::Reset() {
  printf ("Reset Radiosity Solution\n");
  Cleanup();
  Initialize();
}


void Radiosity::Initialize() {
  // create and fill the data structures
  n = mesh->numFaces();
  patches = new Face*[n];
  //formfactors = new float[n*n];
  area = new float[n];
  undistributed = new Vec3f[n];
  absorbed = new Vec3f[n];
  radiance = new Vec3f[n];


  // number the patches, store them in an array (since the bag is unordered)
  Iterator<Face*> *iter = mesh->getFaces()->StartIteration();
  int i = 0;
  while (Face *f = iter->GetNext()) {
    patches[i] = f;
    f->setRadiosityPatchIndex(i);
    setArea(i,patches[i]->getArea());
    i++;
  }
  mesh->getFaces()->EndIteration(iter);
  printf("Just before computing form factors\n");
  //scanf("%*c");
  // compute the form factors and initialize the energies
  for (i = 0; i < n; i++) {
    computeFormFactors(i);
    Vec3f emit = getPatch(i)->getEmit();
    setUndistributed(i,emit);
    setAbsorbed(i,Vec3f(0,0,0));
    setRadiance(i,emit);
  }

  

  // find the patch with the most undistributed energy
  findMaxUndistributed();
  
   
}

void Radiosity::Cleanup() {
  delete [] patches;
  delete [] formfactors;
  delete [] area;
  delete [] undistributed;
  delete [] absorbed;
  delete [] radiance;
}

// ================================================================
// ================================================================

void Radiosity::findMaxUndistributed() {
  // find the patch with the most undistributed energy 
  // don't forget that the patches may have different sizes!
  max_undistributed_patch = -1;
  total_undistributed = 0;
  total_area = 0;
  float max = -1;
  for (int i = 0; i < n; i++) {
    float m = getUndistributed(i).Length() * getArea(i);
    total_undistributed += m;
    total_area += getArea(i);
    if (max < m) {
      max = m;
      max_undistributed_patch = i;
    }
  }
  assert (max_undistributed_patch >= 0 && max_undistributed_patch < n);
}

// =======================================================================
// =======================================================================

bool Radiosity::Iterate() {

  printf ("Iterate\n");

  // compute the radiosity solution!
  
  // return true when the solution error is < epsilon
  // so the animation loop can stop

  return false;
}


void Radiosity::computeFormFactors(int i) {

  // compute the form factors

  // if your method for determining the form factors is only an
  // estimate, you may need to normalize something so you don't
  // gain/lose energy in the system

}

// =======================================================================
// PAINT
// =======================================================================

Vec3f ComputeNormal(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) {
  Vec3f v12 = p2;
  v12 -= p1;
  Vec3f v23 = p3;
  v23 -= p2;
  Vec3f normal;
  Vec3f::Cross3(normal,v12,v23);
  normal.Normalize();
  return normal;
}

Vec3f ComputeNormal(Face *f) {
  return ComputeNormal((*f)[0]->get(),(*f)[1]->get(),(*f)[2]->get());
}

inline float tone_func(float x) {
  assert (x >= 0);
  // a tone mapping hack to map [0,oo) -> [0,1]
  // there's nothing special or physically based about this function
  
  float answer = -exp(-3*x)+1; 
  assert (answer >= 0 && answer <= 1);
  return answer;
}

void Radiosity::insertColor(Vec3f v) {
  if (args->tone_map) {
    float r = tone_func(v.x());
    float g = tone_func(v.y());
    float b = tone_func(v.z());
    glColor3f(r,g,b);
  } else {
    glColor3f(v.x(),v.y(),v.z());
  }
}

void insertNormal(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) {
  Vec3f normal = ComputeNormal(p1,p2,p3);
  glNormal3f(normal.x(), normal.y(), normal.z());
}

Vec3f Radiosity::whichVisualization(enum RENDER_MODE mode, Face *f, int i)
{
  assert (getPatch(i) == f);
  assert (i >= 0 && i < n);
  if (mode == RENDER_LIGHTS) {
    return f->getEmit();
  } else if (mode == RENDER_UNDISTRIBUTED) { 
    return getUndistributed(i);
  } else if (mode == RENDER_ABSORBED) {
    return getAbsorbed(i);
  } else if (mode == RENDER_RADIANCE) {
    return getRadiance(i);
  } else if (mode == RENDER_FORM_FACTORS) {
    float scale = 0.2 * total_area/getArea(i);
    
    float factor = scale * getFormFactor(max_undistributed_patch,i);
    return Vec3f(factor,factor,factor);
  } else {
    assert(0);
  }
}

void CollectFacesWithVertex(Vertex *have, Face *f, Array<Face*> *array) {
  if (array->Member(f)) return;
  if (have != (*f)[0] && have != (*f)[1] && have != (*f)[2]) return;
  array->Add(f);
  for (int i = 0; i < 4; i++) {
    Edge *ea = f->getEdge()->getOpposite();
    Edge *eb = f->getEdge()->getNext()->getOpposite();
    Edge *ec = f->getEdge()->getNext()->getNext()->getOpposite();
   // Edge *ed = f->getEdge()->getNext()->getNext()->getNext()->getOpposite();
    if (ea != NULL) CollectFacesWithVertex(have,ea->getFace(),array);
    if (eb != NULL) CollectFacesWithVertex(have,eb->getFace(),array);
    if (ec != NULL) CollectFacesWithVertex(have,ec->getFace(),array);
   // if (ed != NULL) CollectFacesWithVertex(have,ed->getFace(),array);
  }
}

void Radiosity::insertInterpolatedColor(int index, Face *f, Vertex *v) {
  Array<Face*> faces = Array<Face*>(10);
  CollectFacesWithVertex(v,f,&faces);
  float total = 0;
  Vec3f color = Vec3f(0,0,0);
  Vec3f normal = ComputeNormal(f);
  for (int i = 0; i < faces.Count(); i++) {
    Vec3f normal2 = ComputeNormal(faces[i]);
    if (normal.Dot3(normal2) < 0.9) continue;
    total += faces[i]->getArea();
    color += faces[i]->getArea() * whichVisualization(RENDER_RADIANCE,faces[i],faces[i]->getRadiosityPatchIndex());
  }
  if(total == 0)
  		   for(;;);
  assert (total > 0);
  color /= total;
  insertColor(color);
}

void Radiosity::PaintSelection(ArgParser *args)
{
	  Vec3f center; mesh->getBoundingBox()->getCenter(center);
  float s = 1/mesh->getBoundingBox()->maxDim();
  glScalef(s,s,s);
  glTranslatef(GLCanvas::xtrans-center.x(),GLCanvas::ytrans-center.y(),GLCanvas::ztrans-center.z());
	 // glNewList(GLCanvas::display_list_selection,GL_COMPILE_AND_EXECUTE);
	//glDisable(GL_LIGHTING);
	printf("Came to render_selection mode\n");
	glPushName(1);
	GLuint pname =1;
	for (int i = 0; i < n; i++) {
	  Face *f = patches[i];
	  //printf("%u",pname);
	  glLoadName(pname++);
	 
	 //  Vec3f color = f->getColor();
	  Vec3f a = (*f)[0]->get();
      Vec3f b = (*f)[1]->get();
      Vec3f c = (*f)[2]->get();
	  glBegin(GL_TRIANGLES);
	 // insertNormal(a,b,c); 
	  glVertex3f(a.x(),a.y(),a.z());
      glVertex3f(b.x(),b.y(),b.z());
      glVertex3f(c.x(),c.y(),c.z()); 
	  glEnd();
	  
	  //glPopName();
     }
	//glEndList();
}
void render_stroke_string(void* font, const char* string)
{
	char* p;
	//float width = 0;

	// Center Our Text On The Screen
    //glPushMatrix();
	// Render The Text
	p = (char*) string;
	while (*p != '\0') glutStrokeCharacter(font, *p++);
	//glPopMatrix();
}
void render_string(void* font, const char* string)
{
	char* p = (char*) string;
	while (*p != '\0') glutBitmapCharacter(font, *p++);
}

void Radiosity::Paint(ArgParser *args) {

  // scale it so it fits in the window
  Vec3f center; mesh->getBoundingBox()->getCenter(center);
  float s = 1/mesh->getBoundingBox()->maxDim();
  glScalef(s,s,s);
    glTranslatef(GLCanvas::xtrans-center.x(),GLCanvas::ytrans-center.y(),GLCanvas::ztrans-center.z());
  
  // this offset prevents "z-fighting" bewteen the edges and faces
  // the edges will always win.
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.1,4.0);
//printf("\nrender mode %d\n",args->render_mode);
  if (args->render_mode == RENDER_MATERIALS) {
    // draw the faces with OpenGL lighting, just to understand the geometry
    // (the GL light has nothing to do with the surfaces that emit light!)
//    printf("I am here\n");
    //scanf("%%d"); // CHANGED 
	 	 // glDisable(GL_BLEND);
	//  glEnable(GL_LIGHT1);
	 // glEnable(GL_LIGHT0);
    glBegin (GL_TRIANGLES);
    for (int i = 0; i < n; i++) {
	      Face *f = patches[i];
      Vec3f color = f->getColor();
	  //color = Vec3f(1,1,0);
      insertColor(color);
      Vec3f a = (*f)[0]->get();
      Vec3f b = (*f)[1]->get();
      Vec3f c = (*f)[2]->get();
      //Vec3f d = (*f)[3]->get();
      insertNormal(a,b,c); 
      glVertex3f(a.x(),a.y(),a.z());
      glVertex3f(b.x(),b.y(),b.z());
      glVertex3f(c.x(),c.y(),c.z());
      //glVertex3f(d.x(),d.y(),d.z());
    }
    glEnd();
	//glEnable(GL_BLEND);
  } else if (args->render_mode == RENDER_SELECTION)
  {
	 // glNewList(GLCanvas::display_list_selection,GL_COMPILE_AND_EXECUTE);
	//glDisable(GL_LIGHTING);
	printf("Came to render_selection mode\n");
	glPushName(1);
	GLuint pname =1;
	for (int i = 0; i < n; i++) {
	  Face *f = patches[i];
	  //printf("%u",pname);
	  glLoadName(pname++);
	 
	 //  Vec3f color = f->getColor();
	  Vec3f a = (*f)[0]->get();
      Vec3f b = (*f)[1]->get();
      Vec3f c = (*f)[2]->get();
	  glBegin(GL_TRIANGLES);
	 // insertNormal(a,b,c); 
	  glVertex3f(a.x(),a.y(),a.z());
      glVertex3f(b.x(),b.y(),b.z());
      glVertex3f(c.x(),c.y(),c.z()); 
	  glEnd();
	  
	  //glPopName();
     }
	//glEndList();
	
	//glEnable(GL_LIGHTING);
  } else if (args->render_mode == RENDER_RADIANCE && args->interpolate == true) {
    // interpolate the radiance values with neighboring faces having the same normal
    glDisable(GL_LIGHTING);
    glBegin (GL_TRIANGLES);
    for (int i = 0; i < n; i++) {
      Face *f = patches[i];
      Vec3f a = (*f)[0]->get();
      Vec3f b = (*f)[1]->get();
      Vec3f c = (*f)[2]->get();
     // Vec3f d = (*f)[3]->get();
      insertInterpolatedColor(i,f,(*f)[0]);
      glVertex3f(a.x(),a.y(),a.z());
      insertInterpolatedColor(i,f,(*f)[1]);
      glVertex3f(b.x(),b.y(),b.z());
      insertInterpolatedColor(i,f,(*f)[2]);
      glVertex3f(c.x(),c.y(),c.z());
     // insertInterpolatedColor(i,f,(*f)[3]);
     // glVertex3f(d.x(),d.y(),d.z());
    }
    glEnd();
    glEnable(GL_LIGHTING);
  } else {
		
    // for all other visualizations, just render the patch in a uniform color
//    glDisable(GL_LIGHTING);
//    glBegin (GL_TRIANGLES);
////    printf("%d\n",n);
////    scanf("%*d");
//    for (int i = 0; i < n; i++) {
//      Face *f = patches[i];
//      Vec3f color = whichVisualization(args->render_mode,f,i);
//      insertColor(color);
//      Vec3f a = (*f)[0]->get();
//      Vec3f b = (*f)[1]->get();
//      Vec3f c = (*f)[2]->get();
//     // Vec3f d = (*f)[3]->get();
//     if(f->getColor()==Vec3f(0,0,1))
//     {
//      glVertex3f(a.x(),a.y(),a.z());
//      glVertex3f(b.x(),b.y(),b.z());
//      glVertex3f(c.x(),c.y(),c.z());
//	}
//     // glVertex3f(d.x(),d.y(),d.z());
//    }
//    printf("End\n");
//    glEnd();
//    glEnable(GL_LIGHTING);

  }

  if (args->render_mode == RENDER_FORM_FACTORS) {
    // render the form factors of the patch with the most undistributed light
    glDisable(GL_LIGHTING);
    glColor3f(1,0,0);
    Face *t = getPatch(max_undistributed_patch);
    glLineWidth(3);
    glBegin(GL_LINES);
    Vec3f a = (*t)[0]->get();
    Vec3f b = (*t)[1]->get();
    Vec3f c = (*t)[2]->get();
   // Vec3f d = (*t)[3]->get();
    glVertex3f(a.x(),a.y(),a.z());
    glVertex3f(b.x(),b.y(),b.z());
    glVertex3f(b.x(),b.y(),b.z());
    glVertex3f(c.x(),c.y(),c.z());    
    glVertex3f(c.x(),c.y(),c.z());    
   // glVertex3f(d.x(),d.y(),d.z());    
    //glVertex3f(d.x(),d.y(),d.z());    
    glVertex3f(a.x(),a.y(),a.z());
    glEnd();
    glEnable(GL_LIGHTING);
  }

  glDisable(GL_POLYGON_OFFSET_FILL); 
  HandleGLError(); 
  if(args->show_number)
  {
	  GLboolean bb;
	  glGetBooleanv(GL_LIGHT0,&bb);
	  printf("light0 is %d\n",bb);
	  if(bb==false)
		  glEnable(GL_LIGHT0);
	  vector<trio> vect = getMesh()->undo_operations;
	  glColor3f(1,1,0);
	  GLfloat line_width;
	  glGetFloatv(GL_LINE_WIDTH,&line_width);
	  GLfloat params[2];
	  glGetFloatv(GL_LINE_WIDTH_RANGE,params);
	  glLineWidth(params[0]*0.6+params[1]*0.4);
	  GLdouble model[16],proj[16];
	  GLint view[4];
	  glGetIntegerv(GL_VIEWPORT,view);
	  glGetDoublev(GL_MODELVIEW_MATRIX,model);
	  glGetDoublev(GL_PROJECTION_MATRIX,proj);
	  HandleGLError();
	//  printf("projection matrix:\n");
	  /*
	  for(int counter=0; counter<4; counter++)
	  {
		  for(int counter1=0; counter1<4; counter1++)
		  {
			  printf("%lf ",proj[counter*4+counter1]);
		  }
		  printf("\n");
	  }
	  printf("\n");
	  printf("modelview matrix:\n");
	  for(int counter=0; counter<4; counter++)
	  {
		  for(int counter1=0; counter1<4; counter1++)
		  {
			  printf("%lf ",model[counter*4+counter1]);
		  }
		  printf("\n");
	  }*/

	  GLdouble winx,winy,winz;
	//	printf("__________________________________________________________\n");
      while(!vect.empty())
	  {
		  trio temp = vect.back();
		  Vec3f one,two,three;
		  one = getFaceFromID(temp.a)->getCenter();
		  two = getFaceFromID(temp.b)->getCenter();
		  three = getFaceFromID(temp.c)->getCenter();
		  Vec3f av = (one+two+three)/3;
		  one = getFaceFromID(temp.a)->getNormal();
		  two = getFaceFromID(temp.b)->getNormal();
		  three = getFaceFromID(temp.c)->getNormal();
		  Vec3f av1 = (one+two+three)/3;
		  av1.Normalize();
		  Vec3f av2= av-av1*3;
		  char buff[10];
		  sprintf(buff,"n%d",temp.edit_number);
		  
		  glPushMatrix();

		 // gluLookAt(av.x(),av.y(),av.z(),av2.x(),av2.y(),av2.z(),0,0,1);
		  bool test = gluProject(av.x(),av.y(),av.z(),model,proj,view,&winx,&winy,&winz);
		  if(test == false)
			  printf("gluproject failed\n");
		  //printf("winx winy winz %d %d %d\n",int(winx),int(winy),int(winz));
		  //glMatrixMode(GL_MODELVIEW);
		  //glLoadIdentity();
		  //glTranslatef(0,0,-1);
		  //glRasterPos2f(0.5,0.5);
	//	  printf("av2 %f %f %f\n\n\n",av2.x(),av2.y(),av2.z());
		  
		  //printf("final pos = %lf %lf %lf %lf\n",pos[0],pos[1],pos[2],pos[3]);
		  //printf("Av1 %f %f %f\n",av1.x(),av1.y(),av1.z());
		 
		  Vec3f zvec = Vec3f(0,0,-1);
		  Vec3f zv = Vec3f(0,0,1);
		  Vec3f cross;
		  Vec3f::Cross3(cross,zvec,av1);
		  cross.Normalize();
		  Matrix m= Matrix::MakeAxisRotation(cross,-acos(zvec.Dot3(av1)));
		  
	/*	  for(int counter=0; counter<16; counter++)
		  {
			  printf("%lf ",m.Get(counter/4,counter%4));
			  if(counter%4==3)
				  printf("\n");
		  }*/
		  Vec3f new1 = Vec3f(m.Get(0,0),m.Get(1,0),m.Get(2,0));
		  /*if(av1.Dot3(zv)> 0)
			  new1=Vec3f(0,0,0)-new1;
			  */
		  Vec3f new2;
		  /*
		  if(av1.Dot3(zv)<0)*/
			Vec3f::Cross3(new2,Vec3f(1,0,0),av1);
			/*
		  else
		    Vec3f::Cross3(new2,Vec3f(1,0,0),av1);
			*/
			/*
		  printf("Details :\n AV1");av1.Write();
		  printf("\nCross :");cross.Write();
		  printf("\nnew1 :");new1.Write();
		  printf("\n  new2 : ");new2.Write();
		  printf("\n zv : ");zv.Write();
		  printf("angle of rotation matrix %lf\n",-180/M_PI*acos(zvec.Dot3(av1)));
		 
		  glColor3f(1,0,0);
		  glBegin(GL_LINES);
		  glVertex3f(av.x(),av.y(),av.z());
		  glVertex3f(av.x()+10,av.y(),av.z());
		  glColor3f(0,1,1);
		  glVertex3f(av.x(),av.y(),av.z());
		  glVertex3f(av.x(),av.y(),av.z()+30);
		  glColor3f(1,1,0);
		  glVertex3f(av.x(),av.y(),av.z());
		  glVertex3f(av.x()-10*av1.x(),av.y()-10*av1.y(),av.z()-10*av1.z());
		  glColor3f(1,0,1);
		  glVertex3f(av.x(),av.y(),av.z());
		  glVertex3f(av.x()+10*new2.x(),av.y()+10*new2.y(),av.z()+10*new2.z());
		  glColor3f(1,1,1);
		  glVertex3f(av.x(),av.y(),av.z());
		  glVertex3f(av.x()+100*new1.x(),av.y()+100*new1.y(),av.z()+100*new1.z());
		  glEnd();*/
			new2.Normalize();
		   glTranslatef(av2.x(),av2.y(),av2.z());
          //if(cross.Length()>0.05)
		  {
			glRotatef(180/M_PI*acos(zvec.Dot3(av1)),cross.x(),cross.y(),cross.z());
			//printf("angle of rot needed = %lf\n\n\n\n",(double)(90-180/M_PI*acos(new2.Dot3(new1))));
		  }
		  glRotatef((90-180/M_PI*acos(new2.Dot3(new1))),0,0,1);
		  glScalef(0.03,0.03,0.03);
		//  glRotatef(180.0/M_PI*acos(av1.x()),1,0,0);
		//  glRotatef(180.0/M_PI*acos(av1.y()),0,1,0);
		  render_stroke_string(GLUT_STROKE_ROMAN,buff);
		  glPopMatrix();
		  vect.pop_back();
	  }
	//printf("__________________________________________________________\n");
	  glLineWidth(line_width);
	  if(bb==false)
		  glDisable(GL_LIGHT0);
  }
  if (args->wireframe) {
    mesh->PaintWireframe(); 
  }
  if(args->votes)
  {
	  glDisable(GL_LIGHTING);
	  mesh->PaintVotes(args->vote_number);
	  
	  glEnable(GL_LIGHTING);
  }
//  mesh->PaintCenterLines();
  if(args->points)
  {
	  mesh->PaintConfidence();
	}
  if(args->volume_rendering && args->render_mode !=RENDER_SELECTION)
  {
	  if(!volume_rendering_setup_complete)
	  {
		  volume_rendering_setup_complete=true;
		  mesh->texture_list_init();
		  HandleGLError();
	  }
	// glEnable(GL_TEXTURE_3D);
	 glDisable(GL_LIGHTING);
	 glDisable(GL_DEPTH_TEST);
	 glEnable(GL_TEXTURE_2D);
	 glEnable(GL_BLEND);
	 glColor4f(1,1,1,0.5);
	 //CHANGED!!!!
	 //PFNGLBLENDEQUATIONPROC glBlendEquation = (PFNGLBLENDEQUATIONPROC) wglGetProcAddress("glBlendEquation");
	 //glBlendEquation(GL_MAX);
	 //glBlendFunc(GL_ONE,GL_ONE);
	 //glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	//printf("%s",glGetString(GL_EXTENSIONS));
	 HandleGLError();
	// gluLookAt(0,0,30,0,0,10,0,1,0);
	 glColor3f(1,0,0);
	 //glTranslatef(60,0,-10);
	 //glTranslatef(1.5,1.5,0);
	 glTranslatef(0.5,0.5,0.5);
	 //printf("\nrdepth = %u\n",mesh->rdepth);
	 glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	 for(int counter=mesh->rdepth-1; counter>=0; counter--)
	 {
		 //printf("%u ",mesh->texturexy[counter]);
		 //glColor3f(1.0,0.0,0.0);
		glBindTexture(GL_TEXTURE_2D, mesh->texturexy[counter]);
		glBegin(GL_QUADS);
		glTexCoord2f(0,0);glVertex3f(0,mesh->rlength/SKIP,counter);
		glTexCoord2f(1,0);glVertex3f(mesh->rwidth/SKIP,mesh->rlength/SKIP,counter);
		glTexCoord2f(1,1);glVertex3f(mesh->rwidth/SKIP,0,counter);
		glTexCoord2f(0,1);glVertex3f(0,0,counter);
		glEnd();
	 }
	 /*
	 glColor4f(1,1,1,1);
	 glBegin(GL_TRIANGLES);
	 glVertex3f(-100,-100,0);
	 glVertex3f(-100,100,0);
	 glVertex3f(100,100,100);
	 glEnd();
	 */
	// glTranslatef(-1.5,-1.5,0);

	 glDisable(GL_BLEND);
	 glDisable(GL_TEXTURE_2D);
	 glEnable(GL_DEPTH_TEST);
	 glEnable(GL_LIGHTING);
		//glCallList(mesh->list);
	mesh->getBoundingBox()->Print();

	 // glCallList(GLCanvas::list);
	//	printf("I was called man!\n");
		HandleGLError();
	}
	//draw_transparent_sphere(10);
	/*GLint glmaxnamestackdepth;
	glGetIntegerv(GL_MAX_NAME_STACK_DEPTH,&glmaxnamestackdepth);
	printf("GL_MAX_NAME_STACK_DEPTH %u\n", glmaxnamestackdepth);
	if(glmaxnamestackdepth >=1)
		printf("niceeeeee :D we have enough space in the name depth stack\n");
		*/
}

void draw_transparent_sphere(GLfloat radius)
{
    glDepthMask(GL_FALSE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glColor4f(1.0, 0.0, 0.0, 0.3);
    DrawSphere(radius,3); // file attached
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);

}
Face * Radiosity::getFaceFromID(int id)
{
	//assert(id>=1);
	return patches[id-1];
}
int processHits (GLint hits, GLuint buffer[]);

Vec3f getCenter1(Face * f)
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
void GLCanvas::show_editing(void)
{
	Array<Face*> arr1(1000);
	Array<Face*> arr2(1000);
	//printf("Entering show_editing\n");
	//int n = picked_faces;
	if(picked_faces == 2)
	{
		/*printf("Picked_faces==2");
		glBegin(GL_LINES);
		Vec3f vec = (faces_picked[0]->getCenter())-Vec3f(0,0,0.1);
		glVertex3f(vec.x(),vec.y(),vec.z());
		vec =(faces_picked[1]->getCenter())-Vec3f(0,0,0.1);
		glVertex3f(vec.x(),vec.y(),vec.z());
		glEnd();*/
		faces_picked[1]->setColor(Vec3f(1,0,0));
	}
	else if (picked_faces == 1)
	{
		/*printf("Picked_faces==1");
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		Vec3f v = (faces_picked[0]->getCenter());
		glTranslatef(v.x(),v.y(),v.z());
		DrawSphere(1,3);
		glPopMatrix();*/
		faces_picked[0]->setColor(Vec3f(1,0,0));
	}
	else if(picked_faces == 3)
	{
		faces_picked[0]->setColor(Vec3f(1,1,1));
		faces_picked[1]->setColor(Vec3f(1,1,1));
		trio tri;
		tri.a = num_picked[0];
		tri.b = num_picked[1];
		tri.c = num_picked[2];
		tri.type = curr_type;
		tri.edit_number = radiosity->getMesh()->undo_operations.size();
		strcpy(tri.annotation,curr_annotation);
		radiosity->getMesh()->undo_operations.push_back(tri);
		while(!radiosity->getMesh()->redo_operations.empty())
			radiosity->getMesh()->redo_operations.pop_back();
		for(int counter=0; counter<3; counter++)
		{
			if(faces_picked[counter]==NULL)
				printf("faces_picked[%d] is NULL!\n",counter);
		}
//		printf("Picked_faces==3");
		double h = (faces_picked[0]->getCenter()-faces_picked[1]->getCenter()).Length();
		Vec3f nor = (faces_picked[0]->getCenter())-(faces_picked[1]->getCenter());
		nor.Normalize();
		Vec3f temppoint = ((faces_picked[2]->getCenter())-(faces_picked[1]->getCenter()));
		double lambda = temppoint.Dot3(nor);
		double radius = ((faces_picked[2]->getCenter())-((faces_picked[1]->getCenter())+lambda*nor)).Length();
		
//		printf("About to do BFS\n");
		arr1.Clear();
		arr2.Clear();
		BFS_debug(faces_picked[0],30,&arr1);
		BFS_debug(faces_picked[1],30,&arr2);
//		printf("Reached till the end of doing debug_BFS\n");
		//return;
		for(int counter=0; counter<arr1.Count(); counter++)
		{
			Vec3f center1 = (arr1[counter]->getCenter());
			double lambda_temp = (center1- (faces_picked[1]->getCenter())).Dot3(nor);
			if(lambda_temp >=0 && lambda_temp <= h)
			{
				// its probably inside the cylinder
				double radius_temp = (center1-((faces_picked[1]->getCenter())+lambda_temp*nor)).Length();
				if(radius_temp < radius)
				{
					// its inside
					if(curr_type ==0)
						arr1[counter]->setColor(Vec3f(0,1,0));
					else
						arr1[counter]->setColor(Vec3f(1,0,0));
				}
			}
		}
//		printf("Reached till the end of first loop\n");
		for(int counter=0; counter<arr2.Count(); counter++)
		{
			Vec3f center1 = (arr2[counter]->getCenter());
			double lambda_temp = (center1- (faces_picked[1]->getCenter())).Dot3(nor);
			if(lambda_temp >0 && lambda_temp < h)
			{
				// its probably inside the cylinder
				double radius_temp = (center1-((faces_picked[1]->getCenter())+lambda_temp*nor)).Length();
				if(radius_temp < radius)
				{
					// its inside
					if(curr_type ==0)
						arr2[counter]->setColor(Vec3f(0,1,0));
					else
						arr2[counter]->setColor(Vec3f(1,0,0));
				}
			}
		}
//		printf("Reached till the end of second loop\n");
	}
	//printf("Exiting show_editing\n");
}
void GLCanvas::remove_editing()
{
		Array<Face*> arr1(1000);
		Array<Face*> arr2(1000);
		if(radiosity->getMesh()->undo_operations.empty())
		{
			printf("Nothing to undo!\n");
			return;
		}
		trio top = radiosity->getMesh()->undo_operations.back();
		faces_picked[0]=radiosity->getFaceFromID(top.a);
		faces_picked[1]=radiosity->getFaceFromID(top.b);
		faces_picked[2]=radiosity->getFaceFromID(top.c);
		radiosity->getMesh()->undo_operations.pop_back();
		if(!radiosity->getMesh()->undo_operations.empty())
			strcpy(curr_annotation,radiosity->getMesh()->undo_operations.back().annotation);
		else
		strcpy(curr_annotation,"");
		radiosity->getMesh()->redo_operations.push_back(top);
		double h = (faces_picked[0]->getCenter()-faces_picked[1]->getCenter()).Length();
		Vec3f nor = (faces_picked[0]->getCenter())-(faces_picked[1]->getCenter());
		nor.Normalize();
		Vec3f temppoint = ((faces_picked[2]->getCenter())-(faces_picked[1]->getCenter()));
		double lambda = temppoint.Dot3(nor);
		double radius = ((faces_picked[2]->getCenter())-((faces_picked[1]->getCenter())+lambda*nor)).Length();
		
		//printf("About to do BFS\n");
		arr1.Clear();
		arr2.Clear();
		BFS_debug(faces_picked[0],30,&arr1);
		BFS_debug(faces_picked[1],30,&arr2);
		//printf("Reached till the end of doing debug_BFS\n");
		//return;
		for(int counter=0; counter<arr1.Count(); counter++)
		{
			Vec3f center1 = (arr1[counter]->getCenter());
			double lambda_temp = (center1- (faces_picked[1]->getCenter())).Dot3(nor);
			if(lambda_temp >=0 && lambda_temp <= h)
			{
				// its probably inside the cylinder
				double radius_temp = (center1-((faces_picked[1]->getCenter())+lambda_temp*nor)).Length();
				if(radius_temp < radius)
				{
					// its inside
					arr1[counter]->setColor(Vec3f(1,1,1));
				}
			}
		}
		//printf("Reached till the end of first loop\n");
		for(int counter=0; counter<arr2.Count(); counter++)
		{
			Vec3f center1 = (arr2[counter]->getCenter());
			double lambda_temp = (center1- (faces_picked[1]->getCenter())).Dot3(nor);
			if(lambda_temp >0 && lambda_temp < h)
			{
				// its probably inside the cylinder
				double radius_temp = (center1-((faces_picked[1]->getCenter())+lambda_temp*nor)).Length();
				if(radius_temp < radius)
				{
					// its inside
					arr2[counter]->setColor(Vec3f(1,1,1));
				}
			}
		}
}
void GLCanvas::redo_editing()
{
		Array<Face*> arr1(1000);
		Array<Face*> arr2(1000);
		if(radiosity->getMesh()->redo_operations.empty())
		{
			printf("Nothing to redo!\n");
			return;
		}
		trio top = radiosity->getMesh()->redo_operations.back();
		faces_picked[0]=radiosity->getFaceFromID(top.a);
		faces_picked[1]=radiosity->getFaceFromID(top.b);
		faces_picked[2]=radiosity->getFaceFromID(top.c);
		strcpy(curr_annotation,top.annotation);
		radiosity->getMesh()->redo_operations.pop_back();
		radiosity->getMesh()->undo_operations.push_back(top);

		double h = (faces_picked[0]->getCenter()-faces_picked[1]->getCenter()).Length();
		Vec3f nor = (faces_picked[0]->getCenter())-(faces_picked[1]->getCenter());
		nor.Normalize();
		Vec3f temppoint = ((faces_picked[2]->getCenter())-(faces_picked[1]->getCenter()));
		double lambda = temppoint.Dot3(nor);
		double radius = ((faces_picked[2]->getCenter())-((faces_picked[1]->getCenter())+lambda*nor)).Length();
		
//		printf("About to do BFS\n");
		arr1.Clear();
		arr2.Clear();
		BFS_debug(faces_picked[0],30,&arr1);
		BFS_debug(faces_picked[1],30,&arr2);
//		printf("Reached till the end of doing debug_BFS\n");
		//return;
		for(int counter=0; counter<arr1.Count(); counter++)
		{
			Vec3f center1 = (arr1[counter]->getCenter());
			double lambda_temp = (center1- (faces_picked[1]->getCenter())).Dot3(nor);
			if(lambda_temp >=0 && lambda_temp <= h)
			{
				// its probably inside the cylinder
				double radius_temp = (center1-((faces_picked[1]->getCenter())+lambda_temp*nor)).Length();
				if(radius_temp < radius)
				{
					// its inside
					if(curr_type ==0)
					arr1[counter]->setColor(Vec3f(0,1,0));
					else
					arr1[counter]->setColor(Vec3f(1,0,0));
				}
			}
		}
//		printf("Reached till the end of first loop\n");
		for(int counter=0; counter<arr2.Count(); counter++)
		{
			Vec3f center1 = (arr2[counter]->getCenter());
			double lambda_temp = (center1- (faces_picked[1]->getCenter())).Dot3(nor);
			if(lambda_temp >0 && lambda_temp < h)
			{
				// its probably inside the cylinder
				double radius_temp = (center1-((faces_picked[1]->getCenter())+lambda_temp*nor)).Length();
				if(radius_temp < radius)
				{
					// its inside
					if(curr_type ==0)
					arr2[counter]->setColor(Vec3f(0,1,0));
					else
					arr2[counter]->setColor(Vec3f(1,0,0));
				}
			}
		}
}
void GLCanvas::get_annotation(void)
{
	printf("\nEnter the annotation for the current edit :");
	if( fgets(curr_annotation, sizeof(curr_annotation), stdin) == NULL )
    {
    cerr << "fgets returned null..." << endl;
    }
}
void GLCanvas::load_edits(void)
{
	printf("\n Enter the filename to load the edits from: ");
	char filename[512];
	if( fgets(filename, sizeof(filename), stdin) == NULL )
    {
    cerr << "fgets returned null..." << endl;
    }
	//char ch;
	//while((ch=getc(stdin))!=EOF);
	FILE * fp = fopen(filename,"r");
	if(fp==NULL)
	{
		printf("Couldn't open %s for reading\n",filename);
		return;
	}
	int num;
	trio t;
	int pc = 0;
	while(fscanf(fp,"Edit %d\n",&num)>0)
	{
		pc++;
		if( fscanf(fp,"Type %d\n",&curr_type) == EOF )
      {
      cerr << "EOF encountered by fscanf" << endl;
      }
		t.type = curr_type;
		t.edit_number = num;
		if( fscanf(fp,"Triangles %d %d %d\n",&t.a,&t.b,&t.c) == EOF )
      {
      cerr << "EOF encountered by fscanf" << endl;
      }
		if( fgets(filename,512,fp) == NULL )
      {
      cerr << "fgets returned null..." << endl;
      }
		printf("f-%s\n",filename);
		if( fgets(filename,512,fp) == NULL )
      {
      cerr << "fgets returned null..." << endl;
      }
		printf("f-%s\n",filename);
		if( fgets(filename,512,fp) == NULL )
      {
      cerr << "fgets returned null..." << endl;
      }
		printf("f-%s\n",filename);
		//finished throwing junk
		int unused = fscanf(fp,"Annotation:'");
    unused++;
		if( fgets(t.annotation,512,fp) == NULL )
      {
      cerr << "fgets returned null..." << endl;
      }
	//	printf("\nt.annotation size = %d\n t.annotation = |%s|\n",strlen(t.annotation),t.annotation);
		t.annotation[strlen(t.annotation)-2]=0;
	//	printf("annot-%s-",t.annotation);
		strcpy(curr_annotation,t.annotation);
		Array<Face*> arr1(1000);
		Array<Face*> arr2(1000);

		faces_picked[0]=radiosity->getFaceFromID(t.a);
		faces_picked[1]=radiosity->getFaceFromID(t.b);
		faces_picked[2]=radiosity->getFaceFromID(t.c);
		//radiosity->getMesh()->redo_operations.pop();
		radiosity->getMesh()->undo_operations.push_back(t);
		double h = (faces_picked[0]->getCenter()-faces_picked[1]->getCenter()).Length();
		Vec3f nor = (faces_picked[0]->getCenter())-(faces_picked[1]->getCenter());
		nor.Normalize();
		Vec3f temppoint = ((faces_picked[2]->getCenter())-(faces_picked[1]->getCenter()));
		double lambda = temppoint.Dot3(nor);
		double radius = ((faces_picked[2]->getCenter())-((faces_picked[1]->getCenter())+lambda*nor)).Length();
		
//		printf("About to do BFS\n");
		arr1.Clear();
		arr2.Clear();
		BFS_debug(faces_picked[0],30,&arr1);
		BFS_debug(faces_picked[1],30,&arr2);
//		printf("Reached till the end of doing debug_BFS\n");
		//return;
		for(int counter=0; counter<arr1.Count(); counter++)
		{
			Vec3f center1 = (arr1[counter]->getCenter());
			double lambda_temp = (center1- (faces_picked[1]->getCenter())).Dot3(nor);
			if(lambda_temp >=0 && lambda_temp <= h)
			{
				// its probably inside the cylinder
				double radius_temp = (center1-((faces_picked[1]->getCenter())+lambda_temp*nor)).Length();
				if(radius_temp < radius)
				{
					// its inside
					if(curr_type ==0)
					arr1[counter]->setColor(Vec3f(0,1,0));
					else
					arr1[counter]->setColor(Vec3f(1,0,0));
				}
			}
		}
//		printf("Reached till the end of first loop\n");
		for(int counter=0; counter<arr2.Count(); counter++)
		{
			Vec3f center1 = (arr2[counter]->getCenter());
			double lambda_temp = (center1- (faces_picked[1]->getCenter())).Dot3(nor);
			if(lambda_temp >0 && lambda_temp < h)
			{
				// its probably inside the cylinder
				double radius_temp = (center1-((faces_picked[1]->getCenter())+lambda_temp*nor)).Length();
				if(radius_temp < radius)
				{
					// its inside
					if(curr_type ==0)
					arr2[counter]->setColor(Vec3f(0,1,0));
					else
					arr2[counter]->setColor(Vec3f(1,0,0));
				}
			}
		}
		/*	if(num==1)
			for(;;);*/
		
	}
//	printf("final num %d\n",num);
	printf("Loaded Edits\n");
	fclose(fp);
	
}
void GLCanvas::save_edits(void)
{
	//save the edits in a file
	//ask for the filename
	printf("\nEnter the filename to save the edits: ");
	char filename[512];
	if( fgets(filename, sizeof(filename), stdin) == NULL )
    {
    cerr << "fgets returned null..." << endl;
    }
	printf("Writing to file %s\n",filename);
	FILE * fp = fopen(filename,"w");
	vector<trio> stk = radiosity->getMesh()->undo_operations;
	//int si = stk.size();
	while(!stk.empty())
	{
		trio temp = stk.back();
		Vec3f one,two,three;
		one = radiosity->getFaceFromID(temp.a)->getCenter();
		two = radiosity->getFaceFromID(temp.b)->getCenter();
		three = radiosity->getFaceFromID(temp.c)->getCenter();
		fprintf(fp,"Edit %d\nType %d\nTriangles %d %d %d\npos1 %lf %lf %lf\npos2 %lf %lf %lf\npos3 %lf %lf %lf\nAnnotation:'%s'\n",temp.edit_number,temp.type,temp.a,temp.b,temp.c,one.x(),one.y(),one.z(),two.x(),two.y(),two.z(),three.x(),three.y(),three.z(),temp.annotation);
		stk.pop_back();
	}
	printf("Edits saved\n");
	fclose(fp);
}
void GLCanvas::PickPaint(int x, int y) {

  //assert_stackdump(HandleGLError());
  //meshes->args->render_vis_mode = 1;  
  //glui->sync_live();
  //Mesh *mesh = radiosity->getMesh();
  // select triangle with mouse click!
	GLint viewport[4];
  glSelectBuffer(BUFSIZE,selectBuf);
  (void) glRenderMode (GL_SELECT);
  //assert_stackdump(HandleGLError());
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  //camera->glPlaceCamera();
  glGetIntegerv(GL_VIEWPORT,viewport);
  gluPickMatrix(x,viewport[3]-y-1,5,5,viewport);
  //double ratio = 1; 
  int w = args->width;
  int h = args->height;
 // printf("w = %d h = %d\n",w,h);
  double aspect = double(w)/double(h);
  double asp_angle = camera->getAngle() * 180/M_PI;
  if (aspect > 1) asp_angle /= aspect;
  gluPerspective(asp_angle, aspect, 1, 100.0);
  glMatrixMode(GL_MODELVIEW);
  //camera->glPlaceCamera();  
  glInitNames();
  //render for select here
  //Render::renderForSelect(meshes);
  bool vol_render = args->volume_rendering;
  bool wire_render = args->wireframe;
  args->volume_rendering=false;
  args->wireframe = false;
  //render here
    enum RENDER_MODE curr_rm = args->render_mode;
  args->render_mode = RENDER_SELECTION;
  //printf("I'm about to call display\n");
  // display 
  display_select(0,0);
  // end display
  //printf("I finished calling display!\n");
  //radiosity->Paint(args);
  args->render_mode = curr_rm;
  args->volume_rendering = vol_render;
  args->wireframe = wire_render;
  //assert_stackdump(HandleGLError());
  int hits;
  HandleGLError();
  // restoring the original projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glFlush();
	
  // returning to normal rendering mode
  hits = glRenderMode(GL_RENDER);
 
  // if there are hits process them
  int id = -1;
  if(hits!=0)
	printf("I got a hit! hits  = %d\n",hits);
  else
	  printf("no hits :(\n");
  if (hits != 0)
    id = processHits(hits,selectBuf);
//  assert_stackdump(HandleGLError());
  PICKED = NULL;
  LAST_MOUSE = Vec3f(0,0,0);
  CENTER = Vec3f(0,0,0);
  if (id >= 0) {
	  
//	  printf("\nand I got an id too! %d\n",id);
	  
    Face *e = radiosity->getFaceFromID(id);
	if(picked_faces <3)
	{
		faces_picked[picked_faces]=e;
		num_picked[picked_faces]=id;
		picked_faces ++;
		show_editing();
		if(picked_faces==3)
		{
			strcpy(curr_annotation,"");
			picked_faces=0;
		}
	}
    //assert_stackdump (e != NULL && e->isATriangle());
    PICKED = (Face*)e;
	Array<Face*>  arr(50);
	/*BFS_debug(e,5,&arr);
	for(int counter=0; counter<arr.Count(); counter++)
		arr[counter]->setColor(Vec3f(1,0,1));
	if(e->getColor()==Vec3f(1,0,0))
	{
		printf("Already picked face :(\n");
	}*/
	//e->setColor(Vec3f(1,0,0));
//    Vec3f vert = (*e)[0]->get();
//    Vec3f normal=Vec3f(0,0,0)- PICKED->getNormal();
//    // find the point on the image plane!
//    // do the matrix stuff
//    GLdouble modelview[16], projmatr[16];
//    GLint viewport[4];
//    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
//    glGetDoublev(GL_PROJECTION_MATRIX, projmatr);
//    glGetIntegerv(GL_VIEWPORT, viewport);
//    y = viewport[3]-y-1;   // it's upside down
//
//    GLdouble modelX, modelY, modelZ;
//    gluUnProject(x, y, 0, modelview, projmatr, viewport, &modelX,
//&modelY, &modelZ);
//
//    Vec3f center = Vec3f(modelX,modelY,modelZ);
//
//    gluUnProject(x, y, 1, modelview, projmatr, viewport, &modelX,
//&modelY, &modelZ);
//
//    Vec3f image_point = Vec3f(modelX,modelY,modelZ);
//    Vec3f direction = image_point - center;
//    direction.Normalize();
//
//    CENTER = center;
//
//    double d1 = vert.Dot3(normal);
//    double d2 = center.Dot3(normal);
//    double d3 = (center+direction).Dot3(normal);
//    
//    double depth = (d2-d1) / (d2-d3);
//   
//    LAST_MOUSE = center + depth*direction;
    //double radius = meshes->args->painting_radius;
    //meshes->getTriangleMesh()->Paint(PICKED,LAST_MOUSE, radius,(density-1)/8.0);

  }
  Render();
  //display();
HandleGLError();
  //rerender_sceneCB(0);
 // assert_stackdump(HandleGLError());
}



int processHits (GLint hits, GLuint buffer[]) {
  int i;
  unsigned int j;
  GLuint names, *ptr, minZ;
  GLuint numberOfNames = 0;
  GLuint *ptrNames = 0;
  
  if (hits <= 0) return -1;
  ptr = (GLuint *) buffer;
  minZ = 0xffffffff;
  for (i = 0; i < hits; i++) {	
    names = *ptr;
    ptr++;
//	printf("*ptr = %u minZ = %u\n",*ptr,minZ);
    if (*ptr < minZ) {
      numberOfNames = names;
      minZ = *ptr;
      ptrNames = ptr+2;
    }
    
    ptr += names+2;
  }
  
//  assert_stackdump (numberOfNames == 1);
  //printf("number of names = %d", numberOfNames);
  ptr = ptrNames;
  for (j = 0; j < numberOfNames; j++,ptr++) {
    return *ptr;
  }
//  printf("I came here\n");
  return 0;
}
// =======================================================================
// =======================================================================



