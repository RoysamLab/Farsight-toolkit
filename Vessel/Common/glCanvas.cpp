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

#include "glCanvas.h"

//#include "marchtetra.cpp"

#include "radiosity.h"
#include "mesh.h"
#include "draw_sphere.h"

#ifndef M_PI
#define M_PI (2*acos(double(0)))
#endif
// ========================================================
// static variables of GLCanvas class

// State of the mouse cursor


int GLCanvas::mouseButton;
int GLCanvas::mouseX;
int GLCanvas::mouseY;
bool GLCanvas::dragged;

//GLfloat ztrans=0,xtrans =0,ytrans=0;

int GLCanvas::display_list_index;
int GLCanvas::display_list_selection;
ArgParser* GLCanvas::args;
Camera* GLCanvas::camera;
Radiosity* GLCanvas::radiosity;

// ========================================================


void GLCanvas::InitLight() {
  // Set the last component of the position to 0 to indicate
  // a directional light source

  GLfloat position[4] = { 30,30,100, 1};
  GLfloat diffuse[4] = { 0.75,0.75,0.75,1};
  GLfloat specular[4] = { 0,0,0,1};
  GLfloat ambient[4] = { 0.2, 0.2, 0.2, 1.0 };

  GLfloat zero[4] = {0,0,0,0};
  glLightfv(GL_LIGHT1, GL_POSITION, position);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
  glLightfv(GL_LIGHT1, GL_AMBIENT, zero);
  glEnable(GL_LIGHT1);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  glEnable(GL_COLOR_MATERIAL);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

  GLfloat spec_mat[4] = {1,1,1,1};
  float glexponent = 30;
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &glexponent);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec_mat);

  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  float back_color[] = { 0.0,0.0,0.0,0.0};
  glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, back_color);
  glEnable(GL_LIGHT1);
 // glEnable(GL_TEXTURE_3D);
    xtrans=0;ytrans=0;ztrans=0;
//	display_select(0,0);  
}


void GLCanvas::display(void)
{
	//printf("I'm in the display function\n");
  // Clear the display buffer, set it to the background color
  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set the camera parameters
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  InitLight(); // light will be a headlamp!

  camera->glPlaceCamera();
   // glTranslatef(xtrans,ytrans,-ztrans);  

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  
  glCallList(display_list_index);
  HandleGLError(); 
   
  // Swap the back buffer with the front buffer to display
  // the scene
  glutSwapBuffers(); 
}

void GLCanvas::display_select(int x,int y)
{
	printf("%u",display_list_selection);
	static int first = true;
	//printf("I'm in the display function\n");
  // Clear the display buffer, set it to the background color
  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set the camera parameters
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  InitLight(); // light will be a headlamp!

  camera->glPlaceCamera();
  //  glTranslatef(xtrans,ytrans,-ztrans);  

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  
  if(first)
  {
	  glNewList(display_list_selection,GL_COMPILE_AND_EXECUTE);
	  radiosity->PaintSelection(args);
	  glEndList();
	  first = false;
  }
  else
  {
	  printf("calling display list now\n");
	glCallList(display_list_selection);
  }
  //printf("before gl error");
  HandleGLError(); 
  // printf("After gl error");
  // Swap the back buffer with the front buffer to display
  // the scene
  glutSwapBuffers(); 
}

// ========================================================
// Callback function for window resize
// ========================================================

void GLCanvas::reshape(int w, int h) {
  // Set the OpenGL viewport to fill the entire window
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);

  // Set the camera parameters to reflect the changes
  camera->glInit(w, h);

  args->width = w;
  args->height = h;
}

// ========================================================
// Callback function for mouse click or release
// ========================================================

void GLCanvas::mouse(int button, int state, int x, int y) {
  // Save the current state of the mouse.  This will be
  // used by the 'motion' function
  mouseButton = button;
  mouseX = x;
  mouseY = y;
  	if(args->edit_mode==true&& state == GLUT_UP && button == GLUT_LEFT_BUTTON && dragged==false)
	{
		PickPaint(x,y);
		return;
	}
	if(state==GLUT_DOWN)
		dragged=false;
	//if(state == GL_DOWN)

}

// ========================================================
// Callback function for mouse drag
// ========================================================

void GLCanvas::motion(int x, int y) {
  // Left button = rotation
  // (rotate camera around the up and horizontal vectors)
dragged = true;
  if (mouseButton == GLUT_LEFT_BUTTON) {
    camera->rotateCamera(0.005*(mouseX-x), 0.005*(mouseY-y));
    mouseX = x;
    mouseY = y;
  }
  // Middle button = translation
  // (move camera perpendicular to the direction vector)
  else if (mouseButton == GLUT_MIDDLE_BUTTON) {
    camera->truckCamera((mouseX-x)*0.05, (y-mouseY)*0.05);
    mouseX = x;
    mouseY = y;
  }
  // Right button = zoom
  // (move camera along the direction vector)
  else if (mouseButton == GLUT_RIGHT_BUTTON) {
    camera->dollyCamera((x-mouseX)*0.05);
    mouseX = x;
    mouseY = y;
  }

  // Redraw the scene with the new camera parameters
  glutPostRedisplay();
}

void GLCanvas::passive_motion(int x, int y)
{
	GLdouble modelview[16], projmatr[16];
    GLint viewport[4];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projmatr);
    glGetIntegerv(GL_VIEWPORT, viewport);
    y = viewport[3]-y-1;   // it's upside down
    GLdouble modelX, modelY, modelZ;
    gluUnProject(x, y, 0, modelview, projmatr, viewport, &modelX,&modelY, &modelZ);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslatef(modelX,modelY,modelZ);
	DrawSphere(1,3);
	glPopMatrix();
	glPopMatrix();
}
// ========================================================
// Callback function for keyboard events
// ========================================================

void GLCanvas::keyboard(unsigned char key, int x, int y) {
  switch (key) {
  case 'w':  case 'W':
    args->wireframe = !args->wireframe;
    Render();
    break;
  case 'b':  case 'B':
    args->interpolate = !args->interpolate;
    Render();
    break;
  //case 't':  case 'T':
  //  args->tone_map = !args->tone_map;
  //  Render();
    //break;
  //case 's': case 'S':
  //  radiosity->Cleanup();
  //  radiosity->getMesh()->Subdivision();
  //  radiosity->Initialize();
  //  Render();
  //  break;
  case 'v': case 'V':
    args->render_mode = RENDER_MODE((args->render_mode+1)%NUM_RENDER_MODES);
    switch (args->render_mode) {
    case RENDER_MATERIALS: printf ("RENDER_MATERIALS\n"); break;
    case RENDER_LIGHTS: printf ("RENDER_LIGHTS\n"); break;
    case RENDER_UNDISTRIBUTED: printf ("RENDER_UNDISTRIBUTED\n"); break;
    case RENDER_ABSORBED: printf ("RENDER_ABSORBED\n"); break;
    case RENDER_RADIANCE: printf ("RENDER_RADIANCE\n"); break;
    case RENDER_FORM_FACTORS: printf ("RENDER_FORM_FACTORS\n"); break;
    default: assert(0); }
    Render();
    break;
  //case 'r': case 'R':
  //  radiosity->Reset();
  //  Render();
  //  break;
  case 'o':case 'O':
	  printf("I came here\n");
	  args->votes = !args->votes;
	  Render();
		break;
  case '[': 
		args->vote_number++;
		Render();
			break;
  case ']': 
	  args->vote_number--;
	  Render();
			break;
  case 'p': case 'P':
		args->points = !args->points;
		if(args->points)
		args->render_mode=RENDER_LIGHTS;
		else
		args->render_mode=RENDER_MATERIALS;
		Render();
		break;
  case ' ':
		args->volume_rendering = !args->volume_rendering;
		Render();
		Render();
		break;
  case 'i': case 'I':
    radiosity->Iterate();
    Render();
    break;
  //case 'a': case 'A':
  //  args->animate = !args->animate;
  //  if (args->animate) 
  //    printf ("animation started, press 'A' to stop\n");
  //  else
  //    printf ("animation stopped, press 'A' to start\n");
  //  break;
  case 'q':  case 'Q':
    exit(0);
    break;
	
	/* dont update xtrans ytrans and ztrans.. it causes trouble :-?*/
  case '1':
		 ztrans -= 0.1;
		// printf("going back!\n");
		 break;
  case '2':
  	  ztrans +=0.1;
  	   break;
  case '3':
		 xtrans -= 0.001;
		// printf("going back!\n");
		 break;
  case '4':
  	  xtrans +=0.001;
  		break;
   case '5':
		 ytrans -= 0.001;
		// printf("going back!\n");
		 break;
  case '6':
  	  ytrans +=0.001;
  	   break;
	  
 // case 'e':
	//  args->edit_mode = !args->edit_mode;
	//  if(args->edit_mode == true)
	//  {
	//	  printf("Entering edit mode\n");
	////	  printf("Enabling passive motion function\n");
	////	  glutPassiveMotionFunc(NULL);
	//  }
	//  else
	//  {
	////	  printf("Disabling passive motion function\n");
	////	  glutPassiveMotionFunc(NULL);
	//	  printf("Leaving edit mode\n");

	//  }
	//  break;
  //case 'u': case 'U':
	 // remove_editing();
	 // Render();
	 // break;
  //case 'r': case 'R' :
	 // redo_editing();
	 // Render();
	 // break;
  //case 'a': case 'A':
	 // get_annotation();
	 // break;
  //case 's': case 'S':
	 // save_edits();
	 // break;
  //case 't': case 'T':
	 // curr_type = curr_type+1;
	 // curr_type = curr_type%2;
	 // printf("Toggled edit type to %d\n",curr_type);
	 // break;
  //case 'l':case 'L':
	 // load_edits();
	 // Render();
	 // break;
  //case 'n':case 'N':
	 // args->show_number=!args->show_number;
	 // Render();
	 // break;
  default:
    printf("UNKNOWN KEYBOARD INPUT  '%c'\n", key);
  }
  	glutPostRedisplay();
}

void GLCanvas::idle() {
//  if (args->animate) {
//    if (radiosity->Iterate()) {
//      // when iterate returns true (solution is close enough) stop
//      args->animate = false;
//      printf ("animation stopped, press 'A' to start\n");
//    }
//    Render();
//  }

////////// **************************
//int num = 10;
//	static Bag<Face *>*bf = radiosity->getMesh()->getFaces();
//	static Iterator<Face*>* iter ;
//	static Array<Face*> arr(50);
//	static bool started = false;
//	
//	if(started == false)
//		iter = bf->StartIteration();
//	double maxc = 0.12,minc = -0.06;
//	double median = 0.02;
////	vector <double> curvs;
////	curvs.clear();
//	static int trys = 0;
//	int local_count = 0;
//	Face *f;
//	while(1)
//	{
// 	   f = iter->GetNext();
// 	   if(f == NULL || local_count == 20)
//	    	   break;
//		arr.Clear();
//		BFS(f,num,&arr);
//		//	printf("size %d\n",arr.Count());
//		if(trys%1000==0)
//			printf("%d\n",trys);
//		trys ++;
//		local_count ++;
//		//	if(arr.Count()!=3)
//		//  	printf("something wrong!!!");
//
//		double sum_curvature=0;
//		int sum_count=0;
//		for (int co = 0; co <arr.Count();co++)
//		{
//			if(f!=arr[co])
//			{
//				sum_curvature += FindCurvature(f,arr[co]);
//				sum_count++;
//			}
//		}
//		if(sum_count==0)
//			continue;
//		sum_curvature/=sum_count;
////		curvs.push_back(sum_curvature);
////		maxc=MAX(maxc,sum_curvature);
////		minc=MIN(minc,sum_curvature);
//		hashcurv[(int)f]=sum_curvature;
//		Vec3f color = getColorCode(hashcurv[(int)f],median,minc,maxc);
//		//printf("%lf\t",sum_curvature/sum_count);
//	}
////	sort(curvs.begin(),curvs.end());
////	double median = curvs[int(curvs.size()/2)];
//
//	if(f==NULL)
//			bf->EndIteration(iter);
////	DEBUG("ENDED");
//	started = true;

///////////  ************************
}

// ========================================================
// Initialize all appropriate OpenGL variables, set
// callback functions, and start the main event loop.
// This function will not return but can be terminated
// by calling 'exit(0)'
// ========================================================

void GLCanvas::initialize(ArgParser *_args, Radiosity *_radiosity) {

  args = _args;
  radiosity = _radiosity;
  curr_type = 0;
  camera = new PerspectiveCamera(Vec3f(0,0,2),Vec3f(0,0,-1),Vec3f(0,1,0),5 * M_PI/180.0);

  // Set global lighting parameters
  glEnable(GL_LIGHTING);
  glShadeModel(GL_SMOOTH);

  // Set window parameters
  //printf("have i come here just before glutCreateWindow?66");
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
  glEnable(GL_DEPTH_TEST);
  glutInitWindowSize(args->width, args->height);
  glutInitWindowPosition(100,100);
  glutCreateWindow("OpenGL Viewer");
  //printf("have i come here just after glutCreateWindow?66");
  
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glEnable(GL_NORMALIZE);

  // Ambient light
  Vec3f ambColor = Vec3f(0.2,0.2,0.2); 
  GLfloat ambArr[] = { ambColor.x(), ambColor.y(), ambColor.z(), 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambArr);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  glCullFace(GL_BACK);
  //glEnable(GL_CULL_FACE);

  display_list_index = glGenLists(1);
 display_list_selection = glGenLists(1);
  // Initialize callback functions
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  
  glutIdleFunc(idle);

  Render();

  // Enter the main rendering loop
  glutMainLoop();
}


void GLCanvas::Render() {
  glNewList(display_list_index, GL_COMPILE_AND_EXECUTE);

  // =========================================================
  // put your GL drawing calls inside the display list for efficiency
  radiosity->Paint(args);
  
  // =========================================================
  glEndList();
  glutPostRedisplay();
}

// ========================================================
// ========================================================

int HandleGLError() {
  GLenum error;
  int i = 0;
  while ((error = glGetError()) != GL_NO_ERROR) {
    printf ("GL ERROR(%d):  %s\n", i, gluErrorString(error));
    i++;
  }
  if (i == 0) return 1;
  return 0;
}
