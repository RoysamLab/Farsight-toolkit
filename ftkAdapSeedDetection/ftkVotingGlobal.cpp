// ############################################################################################################################################################################
#include"ftkVotingGlobal.h"

double nftkVot::diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks)/CLOCKS_PER_SEC;
	return diffms;
};

int  nftkVot::round_double( double x ) 
{
 int    ix;
 double dx;

 ix = (int) x;
 dx = (double) x - (double) ix;

 if (dx >= 0.5) return(ix+1);
 else return(ix);
};

int  nftkVot::round_double2( double x ) 
{
 int    ix;
 double dx;

 ix = (int) x;
 dx = (double) x - (double) ix;

 if (dx >= 0.5) return(ix+1);
 else return(ix);
};

void nftkVot::stopProgram( void )
{
	int input;
	std::cout << " -> Wating for U: ";
	std::cin >> input;
}