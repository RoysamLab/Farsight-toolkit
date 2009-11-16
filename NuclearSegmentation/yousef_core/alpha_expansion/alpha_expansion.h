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

#ifndef ALPHA_EXPANSION_H
#define ALPHA_EXPANSION_H

#include <iostream>
#include "GraphCut.h"
#include "GCoptimization.h"
#include "Multi_Color_Graph_Learning_2D.h"

void start_alpha_expansion(float* im, unsigned short* seg_im, float* Dterms, int R, int C, int Z, int K);
void alpha_expansion_2d( float *im, float *sublogImg, unsigned short *subclustImg, int R, int C);

#endif

