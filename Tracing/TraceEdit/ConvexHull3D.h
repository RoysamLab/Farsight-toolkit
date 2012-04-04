#ifndef CONVEXHULL3D_H
#define CONVEXHULL3D_H

#include <algorithm>
#include <cstddef>
#include <math.h>
#include <vector>
#include "Face.h"
#include "TraceBit.h"

// copied by Audrey Cheong

/**
* @author Nicolas Barradeau
* http://en.nicoptere.net
*/
class ConvexHull3D
{
public:
	//inline face function? auto by C++
	ConvexHull3D();

	/**
	 * performs a convexhull in 3D
	 * @param	points the Vector3D cloud
	 * @return a series of indices to create the faces of the hull
	 */
	std::vector<int> getBoundaryPoints(std::vector<TraceBit> &points);
	//static TraceBit getCentroid(std::vector<TraceBit> &points, int index, Face * face);
	void removeDuplicates(std::vector<int> &vec);

private:
	std::vector<Face*> validFaces;
	std::vector<Face*> visibleFaces;
	std::vector<Face*> tempFaces;
};
#endif