 #include "dhSeedGrid.h"

namespace dh
{

SeedGrid::SeedGrid( Histogram * h_in, bool light_bkd, int sample_dist_in )
	: grid_array ( (int)(histSize/(float)sample_dist_in*2+1),
		           (int)(histSize/(float)sample_dist_in*2+1),
				   (int)(histSize/(float)sample_dist_in*2+1) ),
		center( (int)(histSize/(float)sample_dist_in),
				(int)(histSize/(float)sample_dist_in), 
				(int)(histSize/(float)sample_dist_in) )
{ 
	h = h_in;
	sample_dist = sample_dist_in;
	size = (int)( histSize / (float)sample_dist_in*2+1);

	g = grid_array.a;
		
	// Fill Grid
	std::cout << "    Filling grid..." << std::endl;
	int i, j, k, ip, jp, kp, hi, hj, hk;
			
	for (i=0, ip=-(size-1)/2; i<size; i++, ip++)
	{ 
		hi = h->mode[0] + ip*sample_dist;
		if ( hi >= h->d1min && hi <= h->d1max )
		{ 
			for (j=0, jp=-(size-1)/2; j<size; j++, jp++)
			   { 
				hj = h->mode[1] + jp*sample_dist;
			       if ( hj >= h->d3min && hj <= h->d3max )
				   { 
					for (k=0, kp=-(size-1)/2; k<size; k++, kp++)
				       { 
						hk = h->mode[2] + kp*sample_dist;
			               if ( hk >= h->d2min && hk <= h->d2max
							// Point is darker than backgound
							 && (!light_bkd || hk <= h->mode[2]) )
						{ 
							g[i][j][k] = ( h->v(hi, hj, hk) == 0 ) ? 0 : 1;
							//if (g[i][j][k] == 1) DMP( RGB(i, j, k));
						}
					} // end for k
				} // end if hj
			} //end for j
		} // end if hi
	} // end for i
			
	std::cout << "    Finding Surface..." << std::endl;
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
		
void SeedGrid::find_seeds( _RGB& clr1, _RGB& clr2, _RGB& bkgnd )
{ 
	std::cout << "    Finding seeds..." << std::endl;
	if (s.empty())
	{ 
		std::cerr<<"    No surface for seeds!!!!\n";
	}

	double best_val = -1E100;
	_RGB* c1;
	_RGB* c2;
		
	double new_val;
	std::list<_RGB>::iterator it1;

	for( it1 = s.begin(); it1 != s.end(); it1++ )
	{
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
	clr1 = (XYZ)(h->modeAsRGB()) + ((XYZ)(*c1) - (XYZ)center) * sample_dist;
	clr2 = (XYZ)(h->modeAsRGB()) + ((XYZ)(*c2) - (XYZ)center) * sample_dist;
	bkgnd = h->modeAsRGB();
}

//Pass in the seeds and it will change them into most distinct!! 
void SeedGrid::find_most_distinct_colors(_RGB& a1, _RGB& a2, _RGB& bkgrnd)
{
	std::cout << "    Finding most distinct colors..." << std::endl;

	const int max_res = 15;   // MAGIC NUMBER !!! 
	double res[max_res+1];

	// ======= Find Rough Estimate (initial seeds) ===========	
	
	bool a1_moved = false;
	bool a2_moved = false;
	
	do
	{ 
		do
	    { 
			for ( int iter = 0; iter <= max_res; iter++ )
 			{ 
				go_best_dir( a1, a1_moved, a2, bkgrnd );
				go_best_dir( a2, a2_moved, a1, bkgrnd );
				res[iter] = eval_state ( bkgrnd, a1, a2 );
		    }
	    } while ( ( a1_moved || a2_moved ) &&
		           fabs( res[max_res] - res[0] ) >
		           fabs( .02 * res[max_res]) );    // MAGIC NUMBER !!!
      
		int big_jump_size = 8;       // MAGIC NUMBER !!!

		go_best_dir( a1, a1_moved, a2, bkgrnd, big_jump_size, 1 );
		go_best_dir( a2, a2_moved, a1, bkgrnd, big_jump_size, 1 );

	} while ( a1_moved || a2_moved );
}


void SeedGrid::go_best_dir(_RGB& ma, bool& moved, const _RGB& sa, const _RGB& bkgnd, const int r, int res )
{
	double new_val;
	double best_val = NEG_HUGE;
	int best_dR, best_dG, best_dB;
	int x, y, z;
	
	if ( r % res != 0 )
	{ 
		std::cerr<<"    go_best_dir: r % res != 0 may cause instability!" << std::endl; 
	}
	if ( r <= 0 || res <=0 || res > r )
	{ 
		std::cerr<<"    go_best_dir: Invalid parameters: r=" << r << " res=" << res << std::endl;
	}

	for ( x = -r; x <= r; x+= res )
	{ 
		for ( y = -r; y <= r; y+= res )
		{ 
			for ( z = -r; z <= r; z+= res )
			{
				if ( // Point is inside blob                 v Blob Edge Threshold
			        h->v(ma.R+x, ma.G+y, ma.B+z) > 0
					  // and point is darker than background on each axis
				//	  && ma.R+x < at.bkgrnd.R
				//	  && ma.G+y < at.bkgrnd.G 
				//	  && ma.B+z < at.bkgrnd.B    // Intensity, for RLI
				// 4/16/98- this is wrong place for this, with new seeding algorithm
				// add as an option?	 Or make flexible, to work on
				// light on dark OR dark on light?
				   )
			                     
				{ 
					new_val = eval_state ( bkgnd, _RGB(ma.R+x, ma.G+y, ma.B+z), sa );
					if ( new_val > best_val )
					{ 
						best_val = new_val;
						best_dR = x;
						best_dG = y;
						best_dB = z;
					}
				}
			}
		}
	}
	if ( best_val != DBL_MIN )
	{ 
		ma.R += best_dR;
		ma.G += best_dG;
		ma.B += best_dB;
	}
	else
	{ 
		std::cerr<<"    Optimization Algorithm Failure: point is outside blob!!!\n"; 
	}

	moved = ( best_dR == 0 || best_dG == 0 || best_dB == 0 ) ? false : true;
 }

} //end namespace dh
