#ifndef FACE_H
#define FACE_H

#include <math.h>
#include <vector>
#include "TraceBit.h"

// copied by Audrey Cheong

/**
* @author Nicolas Barradeau
* http://en.nicoptere.net
*/

class Face
{
public:
	Face(int index_0, int index_1, int index_2);
	void computePlane();
	bool isVisible(TraceBit const& point);
	TraceBit getCentroid();
	TraceBit getCentroid( int index );
	void flip();
	bool operator ==( Face const* other ) const;
	//bool operator !=( Face const* other ) const;

	static std::vector<TraceBit> points;
	int index_0, index_1, index_2;
	TraceBit point_0, point_1, point_2;
	double a, b, c, d;

private:
	void update();
};
#endif