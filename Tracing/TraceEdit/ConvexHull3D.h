//#ifndef CONVEXHULL3D_H
//#define CONVEXHULL3D_H
//
//#include <math>
//#include <vector>
//#include "TraceBit.h"
//
///**
//* @author Nicolas Barradeau
//* http://en.nicoptere.net
//*/
//
//class ConvexHull3D
//{
//public:
//	//inline face function? auto by C++
//	ConvexHull3D();
//
//	/**
//	 * performs a convexhull in 3D
//	 * @param	points the Vector3D cloud
//	 * @return a series of indices to create the faces of the hull
//	 */
//	std::vector<TraceBit*> process(std::vector<TraceBit*> points);
//	static TraceBit centroid(vector<TraceBit*> points, int index, Face face);
//
//private:
//	std::vector<Face*> validFaces;
//	std::vector<Face*> visibleFaces;
//	std::vector<Face*> tempFaces;
//}
//
//class Face
//{
//public:
//	Face(int i0,int i1,int i2);
//	void computePlane();
//	bool isVisible(TraceBit &p);
//	TraceBit getCentroid();
//	void flip();
//
//	static vector<TraceBit*> points;
//	int i0, i1, i2;
//	int a, b, c, d;
//}
//
//#endif