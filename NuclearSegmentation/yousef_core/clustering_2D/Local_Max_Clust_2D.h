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

#ifndef LOCAL_MAX_CLUST_2D_H
#define LOCAL_MAX_CLUST_2D_H

#include <iostream>
#include <algorithm>

int local_max_clust_2D(float* im_vals, int r, int c, double scale, unsigned short* out1, int* seeds_x, int* seeds_y, int num_seeds, unsigned short* bin_image);

#endif