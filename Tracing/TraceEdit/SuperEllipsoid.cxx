/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#include "SuperEllipsoid.h"

SuperEllipsoid::SuperEllipsoid()
{
}
void SuperEllipsoid::calculateEllipsoid(std::vector<int> &boundary_points)
{
	//TraceBit rootbit = boundary_points[0];
	//double Max = 0;
	////double Min = 1000;
	////int minIndex = -1;
	//int maxIndex = -1;
	//int index = 0;

	//TraceBit components;
	//components.x = 0;
	//components.y = 0;
	//components.z = 0;

	/*
	1) find longest axis by calculated tip to soma euclidean distance //redundant from CellTrace?
	2) get box in the plane formed by soma and long axis,
		a)expand box, increase shape curvature until distance to surface is small
	3) get box at perpendicular axis, apply superellipse
	4) get 3D surface from spherical product of	the 2D 
	*/

	/*
	1) To reduce computation time, work with boundary points only:
		convexhull3D
	2) 
	*/

	/*std::vector<TraceBit*>::iterator iter;
	for (iter = tips.begin(); iter != tips.end(); ++iter, ++index)
	{
		double x = (*iter).x - rootbit.x);
		double y = (*iter).y - rootbit.y);
		double z = (*iter).z - rootbit.z);
		double eucDistance = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
		if (eucDistance > Max)
		{
			Max = eucDistance;
			maxIndex = index;*/
		//}
		//if (eucDistance < Min)
		//{
		//	Min = (*iter);
		//	minIndex = index;
		//}

		////add all vectors
		//components.x += x;
		//components.y += y;
		//components.z += z;
	//}
		

	//double leastSquareDis = 0;

}