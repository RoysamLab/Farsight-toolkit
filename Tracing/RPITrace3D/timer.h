#ifndef TIMER_H
#define TIMER_H

#include <ctime>

class Timer
{
  public:

	Timer() { reset(); }
	void reset()
	{
		start=clock();
	}
	long int elapsedMsec()
	{
		return (clock()-start)*1000/CLOCKS_PER_SEC;
	}

	double elapsedSec()
	{
		return elapsedMsec()/1000 + static_cast<double>(((((clock()-start)*1000/CLOCKS_PER_SEC) % 1000)))/1000.0;
	}

	long int resolutionPerSec()
	{
		return CLOCKS_PER_SEC;
	}

	clock_t start;

};

#endif
