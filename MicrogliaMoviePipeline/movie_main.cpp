#include "SegAndTrace.h"

int main(int argc, char* argv[])
{	
	if( argc == 3)
	{
		MicrogliaMovieSegTracer *segTracer = new MicrogliaMovieSegTracer();
		segTracer->StartSegTracing(argv[1], argv[2]);
	}
	else
	{
		
	}

	return 0;
}