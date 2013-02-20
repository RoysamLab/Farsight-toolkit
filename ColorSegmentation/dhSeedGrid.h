#ifndef _dhSEEDGRID_H_
#define _dhSEEDGRID_H_

#include "dhHistogram.h"
#include "dhEvalState.h"
#include <float.h>
#include <list>

namespace dh
{

class SeedGrid
{ 
public:
	Histogram * h;	// Histogram to be used
		
	Array3D grid_array;	// Array of sampled grid points
	long int*** g;
	
	std::list<_RGB> s;	// List of surface points
		
	int sample_dist;  // Distance between sample points	
						// MAGIC NUMBER !!! ************************ used to be 4
								// Tune for accuracy (hi) vs. Speed (lo)
								// Must use > 1.
	int size;
	_RGB center;

	SeedGrid( Histogram * h_in, bool light_bkd = false, int sample_dist_in = 4 );

	void find_seeds( _RGB& clr1, _RGB& clr2, _RGB& bkgnd );

	void find_most_distinct_colors(_RGB& a1, _RGB& a2, _RGB& bkgrnd);	//Pass in the seeds and it will change them into most distinct!! 

private:
	void go_best_dir(_RGB& ma, bool& moved, const _RGB& sa, const _RGB& bkgnd, const int r = 1, int res = 1 );
	
};
		
} // end namespace dh
#endif
