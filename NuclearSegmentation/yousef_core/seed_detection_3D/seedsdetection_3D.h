#ifndef SEEDSDETECTION_3D_H
#define SEEDSDETECTION_3D_H

int Seeds_Detection_3D( float* IM, float* IM_out, int* IM_bin, int r, int c, int z, double sigma_min, double sigma_max, double scale_xy, double scale_z,int sampl_ratio,int* bImg, int UseDistMap);

#endif
