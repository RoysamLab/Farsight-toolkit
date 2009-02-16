#include "ftl2d_SeedPoint2D.h"
#include <iostream>
void SeedPoint2D::PrintSelf()	{
	std::cout<<"x\ty\tInt\tsc" <<std::endl;
	std::cout<<x <<"\t"<<y<<"\t"<<Intensity<<"\t"<<Scale<<std::endl;
}
