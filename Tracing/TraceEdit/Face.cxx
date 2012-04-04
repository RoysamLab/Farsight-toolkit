
#include "Face.h"

std::vector<TraceBit> Face::points;

Face::Face(int index_0, int index_1, int index_2)
{
	this->index_0 = index_0;
	this->index_1 = index_1;
	this->index_2 = index_2;

	update();
	computePlane();
}

void Face::update()
{
	this->point_0 = points[ index_0 ];
	this->point_1 = points[ index_1 ];
	this->point_2 = points[ index_2 ];
}

void Face::computePlane()
{
	a = point_0.y * (point_1.z - point_2.z) + point_1.y * (point_2.z - point_0.z) + point_2.y * (point_0.z - point_1.z);
	b = point_0.z * (point_1.x - point_2.x) + point_1.z * (point_2.x - point_0.x) + point_2.z * (point_0.x - point_1.x);
	c = point_0.x * (point_1.y - point_2.y) + point_1.x * (point_2.y - point_0.y) + point_2.x * (point_0.y - point_1.y);
	d = -( point_0.x * ( point_1.y * point_2.z - point_2.y * point_1.z ) + point_1.x * (point_2.y * point_0.z - point_0.y * point_2.z) + point_2.x * (point_0.y * point_1.z - point_1.y * point_0.z) );
}

bool Face::isVisible( TraceBit const& point)
{
	return ( a * point.x + b * point.y + c * point.z + d ) > 0;
}

TraceBit Face::getCentroid()
{
	//centroid of 3 points
	TraceBit centroid;
	centroid.x = ( point_0.x + point_1.x + point_2.x ) / 3;
	centroid.y = ( point_0.y + point_1.y + point_2.y ) / 3;
	centroid.z = ( point_0.z + point_1.z + point_2.z ) / 3;

	return centroid;
}

TraceBit Face::getCentroid( int index )
{
	//centroid of 4 points
	TraceBit point = points[ index ];

	TraceBit centroid;
	centroid.x = ( point.x + point_0.x + point_1.x + point_2.x  ) / 4;
	centroid.y = ( point.y + point_0.y + point_1.y + point_2.y  ) / 4;
	centroid.z = ( point.z + point_0.z + point_1.z + point_2.z  ) / 4;
	return centroid;
}

void Face::flip()
{
	int temp = index_0;
	index_0 = index_1;
	index_1 = temp;

	update();
	computePlane();
}

bool Face::operator ==( Face const* other ) const //check if right
{
	if ( this->index_0 != other->index_0 )
		return false;
	if ( this->index_1 != other->index_1 )
		return false;
	if ( this->index_2 != other->index_2 )
		return false;

	return true;
}

//bool Face::operator !=( Face const* other ) const
//{
//	return !(this == other); //check if right
//}