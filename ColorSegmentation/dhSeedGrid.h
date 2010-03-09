#ifndef _dhSEEDGRID_H_
#define _dhSEEDGRID_H_

#include "dhHistogram.h"
#include <list>

namespace dh
{

class SeedGrid
{ 
public:
	const RGBHistogram * h;	// Histogram to be used
		
	Array3D grid_array;	// Array of sampled grid points
	long int*** g;
	
	std::list<_RGB> s;	// List of surface points
		
	const IntensityType sample_dist;  // Distance between sample points	
	const int size;
	const _RGB center;

	SeedGrid( RGBHistogram * h_in, IntensityType sample_dist_in, bool light_bkd = false )
	: h(h_in), sample_dist(sample_dist_in),
	  grid_array ( (int)(128/(float)sample_dist_in*2+1),
		           (int)(128/(float)sample_dist_in*2+1),
				   (int)(128/(float)sample_dist_in*2+1) ),
				  size ( (int)(128/(float)sample_dist_in*2+1) ),
			      center( (int)(128/(float)sample_dist_in),
						  (int)(128/(float)sample_dist_in), 
						  (int)(128/(float)sample_dist_in) )
	{ 
		g = grid_array.a;
		
		// Fill Grid
		std::cout << "Filling grid..." << std::endl;
		int i, j, k, ip, jp, kp, hi, hj, hk;
			
		for (i=0, ip=-(size-1)/2; i<size; i++, ip++)
		{ 
			hi = h->mode.R + ip*sample_dist;
			if ( hi >= h->rmin && hi <= h->rmax )
			{ 
				for (j=0, jp=-(size-1)/2; j<size; j++, jp++)
			    { 
					hj = h->mode.G + jp*sample_dist;
			        if ( hj >= h->gmin && hj <= h->gmax )
				    { 
						for (k=0, kp=-(size-1)/2; k<size; k++, kp++)
				        { 
							hk = h->mode.B + kp*sample_dist;
			                if ( hk >= h->bmin && hk <= h->bmax
								// Point is darker than backgound
								 && (!light_bkd || hk <= h->mode.B) )
							{ 
								g[i][j][k] = ( h->v(hi, hj, hk) == 0 ) ? 0 : 1;
								//if (g[i][j][k] == 1) DMP( RGB(i, j, k));
							}
						} // end for k
					} // end if hj
				} //end for j
			} // end if hi
		} // end for i
			
		std::cout << "Finding Surface..." << std::endl;
		// Find Surface of sampled blob, fill list
		//int num_list_points = 0;
		for (i=1; i<size-1; i++)
		{ 
			for (j=1; j<size-1; j++)
			{ 
				for (k=1; k<size-1; k++)
				{ 
					//std::cout << "i " << i << "  j " << j << "  k " << k << " : " << g[i][j][k] << std::endl;
					if ( g[i][j][k] == 1
					     && (  g[i+1][j][k] == 0
							|| g[i-1][j][k] == 0
							|| g[i][j+1][k] == 0
							|| g[i][j-1][k] == 0
							|| g[i][j][k+1] == 0
							|| g[i][j][k-1] == 0 )
						 )
					{ 
						g[i][j][k] = 2;
						s.push_back( _RGB(i, j, k) );

					} // end if
				} // end for k
			} // end for j
		} // end for i
	}
		
	void find_seeds( _RGB& clr1, _RGB& clr2, double (*eval_state)(const _RGB&, const _RGB&, const _RGB&) )
	{ 
		std::cout << "Finding seeds..." << std::endl;
		if (s.empty())
		{ 
			std::cerr<<"No surface for seeds!!!!\n";
		}

		double best_val = -1E100;
		_RGB* c1;
		_RGB* c2;
			
		double new_val;
		std::list<_RGB>::iterator it1;

		//for( ListElem<RGB>* p1 = s.list; p1 != NULL; p1 = p1->next )
		for( it1 = s.begin(); it1 != s.end(); it1++ )
		{
			//for (ListElem<RGB>* p2 = p1->next; p2 != NULL; p2 = p2->next)
			std::list<_RGB>::iterator it2 = it1;
			it2++;
			for( ; it2 != s.end(); it2++ )
			{ 
				new_val = eval_state( center, *it1, *it2 );
				   // Add "Darker than Background" clause here, if at all
				if ( new_val > best_val )
				{ 
					best_val = new_val;
					c1 = &(*it1);
					c2 = &(*it2);
					//std::cout << new_val << "  " << *c1 << "  " << *c2 << std::endl;
				}
			}
		}
		clr1 = (XYZ)(h->mode) + ((XYZ)(*c1) - (XYZ)center) * sample_dist;
		clr2 = (XYZ)(h->mode) + ((XYZ)(*c2) - (XYZ)center) * sample_dist;
	}
};
		
} // end namespace dh
#endif
