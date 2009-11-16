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

#ifndef ALPHA_EXPANSION_2D_H
#define ALPHA_EXPANSION_2D_H

#include "GCoptimization.h"
#include "GraphCutConstr.cpp"
#include "GraphCutMex.cpp"
#include "Multi_Color_Graph_Learning_2D.h"

int alpha_expansion_2d( int *im, float *sublogImg, int *subclustImg, int R, int C);

#endif
