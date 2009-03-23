#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif 
//#include <GL/glx.h>
#include <stdio.h>
#include <sys/types.h>
//#include <unistd.h>

void Normalize(float v[3]);
void DrawTriangle(float *v1, float *v2, float *v3);
void Subdivide(float *v1, float *v2, float *v3, int depth);
void DrawSphere(float radius, int detail);
