// ====================================================================
// GLCanvas class by Rob Jagnow.
//
// The GLCanvas can be created with a call to the 'initialize' routine,

// ====================================================================

#ifndef _GL_CANVAS_H_
#define _GL_CANVAS_H_
#include <GL/glut.h>
#include <GL/glext.h>
#include <stdlib.h>

#include "argparser.h"
#include "camera.h"

#define CFH(a,b) (((int)a)<<(b*4))


int HandleGLError();
class Radiosity;

// ====================================================================

class GLCanvas {

private:
 
  // State of the mouse cursor
  static int mouseButton;
  static int mouseX;
  static int mouseY;
  static int curr_type;
  // Callback functions for mouse and keyboard events
  static void display(void);
  static void display_select(int,int);
  static void reshape(int w, int h);
  static void mouse(int button, int state, int x, int y);
  static void motion(int x, int y);
  static void passive_motion(int x, int y);
  static void keyboard(unsigned char key, int x, int y);
  static void idle();
  static void PickPaint(int, int);
  static void show_editing(void);
  static void remove_editing(void);
  static void redo_editing(void);
  static void save_edits(void);
  static void get_annotation(void);
  static void load_edits(void);
  
  static bool dragged;
  static ArgParser *args;
  static Camera *camera;
  static Radiosity *radiosity;

  static void InitLight();

  static int display_list_index;
 static int display_list_selection;

public:
 static GLfloat xtrans,ytrans,ztrans;
  // Constructor and destructor
  GLCanvas(void) { }
 ~GLCanvas(void) {  }

  // Set up the canvas and enter the rendering loop
  // Note that this function will not return but can be
  // terminated by calling 'exit(0)'

  void initialize(ArgParser *_args, Radiosity *_radiosity);
 static void Render();

};



// ====================================================================

#endif
