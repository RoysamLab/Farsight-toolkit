#ifndef LABEL2BORDER_H
#define LABEL2BORDER_H

#include "itkLabelBorderImageFilter.h"
#include "itkImageFileWriter.h"

void label2border(int* img,int c, int r, int z);
void label2border2(int* img,int* bordImg, int c, int r, int z); // Added by Yousef on May 29th 2008
void label2border3(int* img,int* segImg, int c, int r, int z);// Added by Yousef on 9-13-2008

#endif

