#include "ConvexHull3D.h"

ConvexHull3D::ConvexHull3D()
{

}

std::vector<int> ConvexHull3D::getBoundaryPoints(std::vector<TraceBit> &points)
{
	if (points.size() < 4)
	{
		std::vector<int> nothing;
		return nothing;
	}
	Face::points = points;

	//calculate first convex tetrahedron

	//creates a face with the first 3 vertices
	Face * face = new Face( 0, 1, 2 );
	//get the centroid of the tetrahedron formed by the first 4 centroid vertices
	TraceBit vertex = face->getCentroid( 3 );

	if ( face->isVisible( vertex ) )
		face->flip();

	Face * face0 = new Face( 3, face->index_0, face->index_1 );
	
	if ( face0->isVisible( vertex ) )
		face0->flip();
	
	Face * face1 = new Face( 3, face->index_1, face->index_2 );
	
	if ( face1->isVisible( vertex ) )
		face1->flip();
	
	Face * face2 = new Face( 3, face->index_2, face->index_0 );
	if ( face2->isVisible( vertex ) )
		face2->flip();

	//store the tetrahedron faces in the valid faces list
	validFaces.push_back(face);
	validFaces.push_back(face0);
	validFaces.push_back(face1);
	validFaces.push_back(face2);


	//so as we have a convex tetrahedron, we can skip the first 4 points
	for ( int point_vector_index = 4; point_vector_index < points.size(); point_vector_index++ )
	{
		//for each available vertices
		vertex.x = points[ point_vector_index ].x;
		vertex.y = points[ point_vector_index ].y;
		vertex.z = points[ point_vector_index ].z;
		
		//checks the point's visibility from all faces
		int visibleFacesCount = 0;
		for ( int validFaces_vector_index = 0; validFaces_vector_index < validFaces.size(); validFaces_vector_index++ )
		{	
			if ( face->isVisible( vertex ) )
			{
				visibleFaces.push_back( validFaces[validFaces_vector_index] );
				visibleFacesCount++;
			}
		}
		
		//the vertex is not visible : it is inside the convex hull, keep on
		if ( visibleFacesCount == 0 )
		{
			continue;
		}
		
		//the vertex is outside the convex hull
		//delete all visible faces from the valid List
		for ( int visibleFaces_index = 0; visibleFaces_index < visibleFaces.size(); visibleFaces_index++ )
		{
			for ( int validFaces_index = 0; validFaces_index < validFaces.size(); validFaces_index++ )
			{
				if ( visibleFaces[visibleFaces_index] == validFaces[validFaces_index] )
				{
					validFaces.erase(validFaces.begin()+validFaces_index);
					//validFaces.splice( validFaces.indexOf( face ), 1 );
				}
			}
		}
		
		//special case : only one face is visible
		//it's ok to create 3 faces directly for they won't enclose any other point
		if ( visibleFacesCount == 1 )
		{
			face = visibleFaces[ 0 ];
			validFaces.push_back( new Face( point_vector_index, face->index_0, face->index_1 ) );
			validFaces.push_back( new Face( point_vector_index, face->index_1, face->index_2 ) );
			validFaces.push_back( new Face( point_vector_index, face->index_2, face->index_0 ) );
			continue;
		}
		
		//creates all possible new faces from the visibleFaces
		int tempFacesCount = 0;
		for ( int visibleFaces_vector_index = 0; visibleFaces_vector_index < visibleFaces.size(); visibleFaces_vector_index++ )
		{
			tempFaces.push_back( new Face( point_vector_index, visibleFaces[visibleFaces_vector_index]->index_0, visibleFaces[visibleFaces_vector_index]->index_1 ) );
			tempFaces.push_back( new Face( point_vector_index, visibleFaces[visibleFaces_vector_index]->index_1, visibleFaces[visibleFaces_vector_index]->index_2 ) );
			tempFaces.push_back( new Face( point_vector_index, visibleFaces[visibleFaces_vector_index]->index_2, visibleFaces[visibleFaces_vector_index]->index_0 ) );
		}
		
		//test each face against another face in tempFaces
		bool validCurrentIndex;
		for ( int current_index = 0; current_index < tempFaces.size(); current_index++ )
		{
			//search if there is a point in front of the face : 
			//this means the face doesn't belong to the convex hull
			validCurrentIndex = true;
			for ( int other_index = 0; other_index < tempFaces.size(); other_index++)
			{
				if ( current_index != other_index )
				{
					if ( tempFaces[current_index]->isVisible( tempFaces[other_index]->getCentroid() ) )
					{
						validCurrentIndex = false;
						break;
					}
				}
			}
			//the face has no point in front of it
			if ( validCurrentIndex )
			{
				//visibleFaces.push_back( tempFaces[current_index] );
				validFaces.push_back( tempFaces[current_index] ); //crash here, why?
			}
		}
	}
	
	//return indices of the convex hull boundary points
	std::vector<int> result;
	for (int index = 0; index < validFaces.size(); index++)
	{
		result.push_back( validFaces[index]->index_0 );
		result.push_back( validFaces[index]->index_1 );
		result.push_back( validFaces[index]->index_2 );
	}
	removeDuplicates( result );
	return result;
}

void ConvexHull3D::removeDuplicates(std::vector<int> &vec)
{
	//remove duplicate indices
	std::sort( vec.begin(),vec.end() );
	vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );
}