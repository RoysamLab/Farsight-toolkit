/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef SEEDSDETECTION_3D_H
#define SEEDSDETECTION_3D_H

// #include "itkMultiThreader.h"

int Seeds_Detection_3D( float* IM, float** IM_out, unsigned short** IM_bin, size_t r, size_t c, size_t z, double *sigma_min, double *sigma_max, double *scale_xy, double *scale_z,int sampl_ratio, unsigned short* bImg, int UseDistMap, int* minIMout, bool paramEstimation);

#endif
