#ifndef LOCAL_MAX_CLUST_2D_H
#define LOCAL_MAX_CLUST_2D_H

#include <iostream>
#include <algorithm>

int local_max_clust_2D(float* im_vals, int r, int c, double scale, int* out1, int* seeds_x, int* seeds_y, int num_seeds, int* bin_image);

#endif