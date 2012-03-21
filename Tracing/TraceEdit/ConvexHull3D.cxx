#include "ConvexHull3D.h"
//
//ConvexHull3D::ConvexHull3D()
//{
//}
//
//std::vector<int> ConvexHull3D::process(std::vector<TraceBit*> points)
//{
//	if (points.size < 3)
//	{
//		return -1;
//	}
//	Face.points = points;
//
//	//calculate first convex tetrahedron
//			
//	//creates a face with the first 3 vertices
//	Face face = new Face( 0, 1, 2 );
//	//get the centroid of the tetrahedron formed by the first 4 vertices
//	TraceBit v = centroid( points, 3, face);
//
//	if ( face.isVisible( v ) )
//		face.flip();
//
//	Face face0 = new Face( 3, face.i0, face.i1 );
//	
//	if ( face0.isVisible( v ) )
//		face0.flip();
//	
//	Face face1 = new Face( 3, face.i1, face.i2 );
//	
//	if ( face1.isVisible( v ) )
//		face1.flip();
//	
//	Face face2 = new Face( 3, face.i2, face.i0 );
//	if ( face2.isVisible( v ) )
//		face2.flip();
//
//	//store the tetrahedron faces in the valid faces list
//	validFaces = std::vector<Face*>( [ face, face0, face1, face2 ] );
//
//	visibleFaces = new std::vector<Face*>();
//	tempFaces = new std::vector<Face*>();
//
//
//
//	//so as we have a convex tetrahedron, we can skip the first 4 points
//	for ( int i = 4; i < points.size; i++ )
//	{
//		//for each available vertices
//		v = points[ i ];
//		
//		//checks the point's visibility from all faces
//		visibleFaces.size = 0;
//		for ( int j = 0; j < validFaces.size; j++ ) //for each( face in validFaces )
//		{	
//			if ( face.isVisible( v ) )
//			{
//				visibleFaces.push_back( face ); 
//			}
//		}
//		
//		//the vertex is not visible : it is inside the convex hull, keep on
//		if ( visibleFaces.size == 0 )
//		{
//			continue;
//		}
//		
//		//the vertex is outside the convex hull
//		//delete all visible faces from the valid List
//		for ( int j = 0; j < visibleFaces.size; j++ ) //for each ( face in visibleFaces )
//		{
//			validFaces.splice( validFaces.indexOf( face ), 1 );
//		}
//		
//		//special case : only one face is visible
//		//it's ok to create 3 faces directly for they won't enclose any other point
//		if ( visibleFaces.size == 1 )
//		{
//			face = visibleFaces[ 0 ];
//			validFaces.push_back( new Face( i, face.i0, face.i1 ) );
//			validFaces.push_back( new Face( i, face.i1, face.i2 ) );
//			validFaces.push_back( new Face( i, face.i2, face.i0 ) );
//			continue;
//		}
//		
//		//creates all possible new faces from the visibleFaces
//		tempFaces.size = 0;
//		for ( int j = 0; j < visibleFaces.size; j++ ) //for each( face in visibleFaces )
//		{
//			tempFaces.push_back( new Face( i, face.i0, face.i1 ) );
//			tempFaces.push_back( new Face( i, face.i1, face.i2 ) );
//			tempFaces.push_back( new Face( i, face.i2, face.i0 ) );
//		}
//		
//		Face other;
//		for ( int j = 0; j < tmpFaces.size; j++ ) //for each( face in tmpFaces )
//		{
//			//search if there is a point in front of the face : 
//			//this means the face doesn't belong to the convex hull
//			for ( int k = 0; k < tempFaces.size; k++) // search : for each( other in tmpFaces )
//			{
//				if ( face != other )
//				{
//					if ( face.isVisible( other.centroid ) )
//					{
//						face = null;
//						break search;
//					}
//				}
//			}
//			//the face has no point in front of it
//			if ( face != null ) 
//				validFaces.push( face );
//		}
//	}
//	
//	std::vector<int> result;	//var result:Vector.<int> = new Vector.<int>();
//	for (int i = 0; i < validFaces.size; i++)  //for each( face in validFaces ) 
//	{
//		result.push_back( face.i0 );
//
//		result.push( face.i0, face.i1, face.i2 );
//	}
//	return result;
//}
//
//TraceBit ConvexHull3D::centroid(vector<TraceBit*> points, int index, Face face)
//{
//	TraceBit p = points[index];
//	TraceBit p0 = points[ face.i0 ];
//	TraceBit p1 = points[ face.i1 ];
//	TraceBit p2 = points[ face.i2 ];
//
//	TraceBit center;
//	center.x = ( p.x + p0.x +  p1.x +  p2.x  ) / 4;
//	center.y = ( p.y + p0.y +  p1.y +  p2.y  ) / 4;
//	center.z = ( p.z + p0.z +  p1.z +  p2.z  ) / 4;
//	return center;
//}
//
//Face::Face(int i0,int i1,int i2)
//{
//	this->i0 = i0;
//	this->i1 = i1;
//	this->i2 = i2;
//
//	computePlane();
//}
//void Face::computePlane()
//{
//	TraceBit v1 = points[ i0 ];
//	TraceBit v2 = points[ i1 ];
//	TraceBit v3 = points[ i2 ];
//
//	a = v1.y * (v2.z - v3.z) + v2.y * (v3.z - v1.z) + v3.y * (v1.z - v2.z);
//	b = v1.z * (v2.x - v3.x) + v2.z * (v3.x - v1.x) + v3.z * (v1.x - v2.x);
//	c = v1.x * (v2.y - v3.y) + v2.x * (v3.y - v1.y) + v3.x * (v1.y - v2.y);
//	d = -( v1.x * ( v2.y * v3.z - v3.y * v2.z ) + v2.x * (v3.y * v1.z - v1.y * v3.z) + v3.x * (v1.y * v2.z - v2.y * v1.z) );
//}
//bool Face::isVisible(TraceBit &p)
//{
//	return ( a * p.x + b * p.y + c * p.z + d ) > 0;
//}
//TraceBit Face::getCentroid()
//{
//	TraceBit p0 = points[ i0 ];
//	TraceBit p1 = points[ i1 ];
//	TraceBit p2 = points[ i2 ];
//
//	TraceBit centroid;
//	centroid.x = ( p0.x + p1.x + p2.x ) / 3;
//	centroid.y = ( p0.y + p1.y + p2.y ) / 3;
//	centroid.z = ( p0.z + p1.z + p2.z ) / 3;
//
//	return centroid;
//}
//void Face::flip()
//{
//	int temp = i0;
//	i0 = i1;
//	i1 = temp;
//	computePlane();
//}