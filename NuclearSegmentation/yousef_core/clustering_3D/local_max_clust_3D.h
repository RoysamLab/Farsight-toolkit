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

#ifndef LOCAL_MAX_CLUST_3D_H
#define LOCAL_MAX_CLUST_3D_H

#include <iostream>
#include <algorithm>

void local_max_clust_3D(float* im_vals, unsigned short* local_max_vals, unsigned short* bImg, unsigned short* out1, int r, int c, int z, int scale_xy, int scale_z);

#endif
