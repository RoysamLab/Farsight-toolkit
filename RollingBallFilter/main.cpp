#include "RollingBallFilter.h"

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " input_image output_image mean_filter_radius" << std::endl;
		return -1;
	}

	RollingBallFilter* rbf_obj = new RollingBallFilter(argv[1], atof(argv[2]));
	rbf_obj->RunFilter();
	rbf_obj->WriteOutput(argv[3]);
}