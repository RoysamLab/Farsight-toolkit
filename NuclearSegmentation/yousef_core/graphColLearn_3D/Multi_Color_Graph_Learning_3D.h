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

#ifndef Multi_Color_Graph_Learning_3D_H
#define Multi_Color_Graph_Learning_3D_H

#include <vector>
#include <iostream>
#include <math.h>

float* multiColGraphLearning(float* X_vals, unsigned short* labs_vals, unsigned short* color_im, size_t r, size_t c, size_t z, int *NC, int refinemetRange);
//added by Yousef on 11/3/2008
void distToEdge(int *** edge_im, int R, int C, int Z);
//////////////////////////////

#endif
