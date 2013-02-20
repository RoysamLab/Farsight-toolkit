#ifndef _kNearestObjects_txx
#define _kNearestObjects_txx

#include "kNearestObjects.h"

template <int num_dimensions>
kNearestObjects<num_dimensions>::kNearestObjects(std::map< unsigned int, std::vector<double> > centroidMap)
{
	allIds = 0;
	featureTable = NULL;
	this->centerMap = centroidMap;
	sample = SampleType::New();
	treeGenerator = TreeGeneratorType::New();

	sample->SetMeasurementVectorSize( num_dimensions );

	// store the object centroids as measurement vectors and form a map of object ID to centroid
	// also store the measurement vectors as a list sample for the generation of the KD tree

	for (It = centerMap.begin(); It != centerMap.end(); ++It )
	{
		MeasurementVectorType mv;
		unsigned int id;
		id = (*It).first;
		for(int i=0; i<num_dimensions; ++i)
		{
			mv[i] = centerMap[id].at(i);
		}
		//mv[0] = centerMap[id].at(0);
		//mv[1] = centerMap[id].at(1);
		//mv[2] = centerMap[id].at(2);
		idToCentroidMap[id] = mv;
		sample->PushBack( mv );
    }

	
	treeGenerator->SetSample( sample );
	treeGenerator->SetBucketSize( 16 );
	treeGenerator->Update();

	tree = treeGenerator->GetOutput();	
}





//***********************************************************************************************
//***********************************************************************************************
// K Nearest Neighbors
//***********************************************************************************************
//***********************************************************************************************

// returns the k nearest neighbors for all IDs and returns a vector of vectors where
// each inner vector contains an ID itself and its neighbors paired with their distances to the ID
template <int num_dimensions>
std::vector< std::vector< std::pair<unsigned int, double> > > kNearestObjects<num_dimensions>::k_nearest_neighbors_All(unsigned int k, unsigned short Class_dest, unsigned short Class_src)
{	
	// Initialize the vector of vectors
	std::vector<std::vector< std::pair<unsigned int, double> > > kNearestIDs;
	
	// fetch each ID to get its nearest neighbors
	for(IdIt = idToCentroidMap.begin(); IdIt != idToCentroidMap.end(); ++IdIt )
	{
		// initialize a vector to store the neighbors and their distances as pairs
		std::vector< std::pair<unsigned int, double> > kNearestIds;
		unsigned int ID = (*IdIt).first;
		
		// for any ID
		if(Class_src == 0)
		{
			kNearestIds = k_nearest_neighbors_ID(ID, k, Class_dest);
			kNearestIDs.push_back(kNearestIds);
		}
		
		// for ID of a given class
		else
		{
			for(int col=((int)featureTable->GetNumberOfColumns())-1; col>=0; --col)
			{	
				std::string current_column = featureTable->GetColumnName(col);
				if(current_column.find("prediction") != std::string::npos )
				{
					for(int row = 0; row<(int)featureTable->GetNumberOfRows(); ++row)
					{
						if( (ID == featureTable->GetValue(row, 0).ToUnsignedInt())
							&& (featureTable->GetValue(row, col).ToUnsignedShort() == Class_src) )
						{
							kNearestIds = k_nearest_neighbors_ID(ID, k, Class_dest);
							kNearestIDs.push_back(kNearestIds);
							break;
						}
					}
					break;
				}
			}
		}		
	}

	return kNearestIDs;
}



// returns the k nearest neighbors for a set of IDs and returns a vector of vectors where
// each inner vector contains an ID itself and its neighbors paired with their distances to the ID
template <int num_dimensions>
std::vector< std::vector< std::pair<unsigned int, double> > > kNearestObjects<num_dimensions>::k_nearest_neighbors_IDs(std::vector<unsigned int> IDs, unsigned int k, unsigned short Class_dest)
{
	std::vector<std::vector< std::pair<unsigned int, double> > > kNearestIDs;
	for(int n=0; n<(int)IDs.size(); ++n)
	{
		std::vector< std::pair<unsigned int, double> > kNearestIds;
		unsigned int ID = IDs.at(n);
		kNearestIds = k_nearest_neighbors_ID(ID, k, Class_dest);
		kNearestIDs.push_back(kNearestIds);				
	}

	return kNearestIDs;
}



// returns the k nearest neighbors for an ID and returns a vector that contains the
// ID itself and its neighbors paired with their distances to the ID
template <int num_dimensions>
std::vector< std::pair<unsigned int, double> > kNearestObjects<num_dimensions>::k_nearest_neighbors_ID(unsigned int id, unsigned int k, unsigned short Class_dest)
{
	std::vector< std::pair<unsigned int, double> > kNearestIds;
	
	// to get the K nearest neighbors of all classes
	if(Class_dest == 0)
	{
		typename DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
		typename DistanceMetricType::OriginType origin( num_dimensions );
		for (unsigned int i = 0 ; i < (unsigned int)sample->GetMeasurementVectorSize() ; ++i)
		{
			origin[i] = idToCentroidMap[id][i];
		}
		distanceMetric->SetOrigin( origin );

		typename TreeType::InstanceIdentifierVectorType kNeighbors;
		tree->Search( idToCentroidMap[id], k+1, kNeighbors ); 
		kNearestIds.push_back( std::make_pair(id,0) );
		for ( unsigned int i = 0 ; i < kNeighbors.size() ; ++i )
		{
		    typename std::map<int, MeasurementVectorType>::iterator iter;
			MeasurementVectorType mvt = tree->GetMeasurementVector(kNeighbors[i]);
			for(iter = idToCentroidMap.begin(); iter != idToCentroidMap.end(); ++iter)
			{
				if(((*iter).second == mvt) && ((*iter).first!=id))
				{
					double dist = distanceMetric->Evaluate( tree->GetMeasurementVector( kNeighbors[i] ));
					std::pair<unsigned int, double> Neighbor = std::make_pair((*iter).first, dist);
					kNearestIds.push_back(Neighbor);
					break;				
				}
			}
		}
	}

	// to get the K nearest neighbors of a particular class
	else
	{
		typename DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
		typename DistanceMetricType::OriginType origin( num_dimensions );
		for (unsigned int i = 0 ; i < (unsigned int)sample->GetMeasurementVectorSize() ; ++i)
		{
			origin[i] = idToCentroidMap[id][i];
		}
		distanceMetric->SetOrigin( origin );

		kNearestIds.push_back( std::make_pair(id,0) );
	
		int count = 0;
		unsigned int num_neighbors = k+1;
		int increment = 10;
		typename TreeType::InstanceIdentifierVectorType temp_vector;
		
		// repeat till the number of nearest neighbors is greater than or equal to k
		while(count < k)
		{
			count = 0;
			temp_vector.clear();
			typename TreeType::InstanceIdentifierVectorType kNeighbors;
			tree->Search( idToCentroidMap[id], num_neighbors, kNeighbors ); 
			for ( unsigned int i = 0 ; i < kNeighbors.size() ; ++i )
			{
			    typename std::map<int, MeasurementVectorType>::iterator iter;
				MeasurementVectorType mvt = tree->GetMeasurementVector(kNeighbors[i]);
				for(iter = idToCentroidMap.begin(); iter != idToCentroidMap.end(); ++iter)
				{
					// find the ID of the neighbor from the centroid map and check if it is the source ID itself
					if(((*iter).second == mvt) && ((*iter).first!=id))
					{
						for(int col=((int)featureTable->GetNumberOfColumns())-1; col>=0; --col)
						{	
							std::string current_column = featureTable->GetColumnName(col);
							if(current_column.find("prediction") != std::string::npos )
							{
								for(int row = 0; row<(int)featureTable->GetNumberOfRows(); ++row)
								{
									// check if the neighbor is of the desired class
									if( ((*iter).first == featureTable->GetValue(row, 0).ToUnsignedInt())
										&& (featureTable->GetValue(row, col).ToUnsignedShort() == Class_dest) )
									{
										// push the neighbor into a temporary vector and increase the count
										temp_vector.push_back(kNeighbors[i]);
										count++;
										break;
									}
								}
								break;
							}
						}					
						break;				
					}
				}
			}

			increment = increment * 2;;
			num_neighbors += increment;
		}

		// calculate the distance of each neighbor and store it in a vector of pairs
		for ( unsigned int i = 0 ; i < temp_vector.size() ; ++i )
		{
			typename std::map<int, MeasurementVectorType>::iterator iter;
			MeasurementVectorType mvt = tree->GetMeasurementVector(temp_vector[i]);
			{
				for(iter = idToCentroidMap.begin(); iter != idToCentroidMap.end(); ++iter)
				{
					if((*iter).second == mvt)
					{
						double dist = distanceMetric->Evaluate(mvt);
						std::pair<unsigned int, double> Neighbor = std::make_pair((*iter).first, dist);
						kNearestIds.push_back(Neighbor);
						break;
					}
				}
			}
		}
		
		//get the five nearest neighbors from the above vector 
		double max_dist;
		unsigned int max_pos;
		while((unsigned int)kNearestIds.size() != k+1)
		{
			max_dist = 0;
			for(unsigned int x=1; x<k+1; ++x)
			{
				if(kNearestIds.at(x).second > max_dist)
				{
					max_dist = kNearestIds.at(x).second;
					max_pos = x;
				}
			}
			if(kNearestIds.at(k+1).second < max_dist )
				kNearestIds.erase(kNearestIds.begin()+max_pos);
			else
				kNearestIds.erase(kNearestIds.begin()+k+1);
		}
	}

	return kNearestIds;

}





//***********************************************************************************************
//***********************************************************************************************
// Neighbors with certain radius
//***********************************************************************************************
//***********************************************************************************************

// returns the neighbors within radius for all IDs and returns a vector of vectors where
// each inner vector contains an ID itself and its neighbors paired with their distances to the ID
template <int num_dimensions>
std::vector< std::vector< std::pair<unsigned int, double> > > kNearestObjects<num_dimensions>::neighborsWithinRadius_All(double radius, unsigned short Class_dest, unsigned short Class_src)
{	
	// Initialize the vector of vectors
	std::vector<std::vector< std::pair<unsigned int, double> > > inRadiusIDs;

	// fetch each ID to get its nearest neighbors
	for(IdIt = idToCentroidMap.begin(); IdIt != idToCentroidMap.end(); ++IdIt )
	{
		// initialize a vector to store the neighbors and their distances as pairs
		std::vector< std::pair<unsigned int, double> > inRadiusIds;
		unsigned int ID = (*IdIt).first;
		
		// for any ID
		if(Class_src == 0)
		{
			inRadiusIds = neighborsWithinRadius_ID(ID, radius, Class_dest);
			inRadiusIDs.push_back(inRadiusIds);
		}
		
		// for ID of a given class
		else
		{
			for(int col=((int)featureTable->GetNumberOfColumns())-1; col>=0; --col)
			{	
				std::string current_column = featureTable->GetColumnName(col);
				if(current_column.find("prediction") != std::string::npos )
				{
					for(int row = 0; row<(int)featureTable->GetNumberOfRows(); ++row)
					{
						if( (ID == featureTable->GetValue(row, 0).ToUnsignedInt())
							&& (featureTable->GetValue(row, col).ToUnsignedShort() == Class_src) )
						{
							inRadiusIds = neighborsWithinRadius_ID(ID, radius, Class_dest);
							inRadiusIDs.push_back(inRadiusIds);
							break;
						}
					}
					break;
				}
			}
		}		
	}

	return inRadiusIDs;
}



// returns the neighbors within radius for a set of IDs and returns a vector of vectors where
// each inner vector contains an ID itself and its neighbors paired with their distances to the ID
template <int num_dimensions>
std::vector< std::vector< std::pair<unsigned int, double> > > kNearestObjects<num_dimensions>::neighborsWithinRadius_IDs(std::vector<unsigned int> IDs, double radius, unsigned short Class_dest)
{
	std::vector<std::vector< std::pair<unsigned int, double> > > inRadiusIDs;
	for(int n=0; n<(int)IDs.size(); ++n)
	{
		std::vector< std::pair<unsigned int, double> > inRadiusIds;
		unsigned int ID = IDs.at(n);
		inRadiusIds = neighborsWithinRadius_ID(ID, radius, Class_dest);
		inRadiusIDs.push_back(inRadiusIds);				
	}

	return inRadiusIDs;
}



// returns the neighbors within radius for an ID and returns a vector that contains the
// ID itself and its neighbors paired with their distances to the ID
template <int num_dimensions>
std::vector< std::pair<unsigned int, double> > kNearestObjects<num_dimensions>::neighborsWithinRadius_ID(unsigned int id, double radius, unsigned short Class_dest)
{
	typename DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
	typename DistanceMetricType::OriginType origin( num_dimensions );
	for ( unsigned int i = 0 ; i < (unsigned int)sample->GetMeasurementVectorSize() ; ++i )
	{
		origin[i] = idToCentroidMap[id][i];
    }
	distanceMetric->SetOrigin( origin );

	typename TreeType::InstanceIdentifierVectorType radNeighbors;
	tree->Search( idToCentroidMap[id], radius, radNeighbors ); 
	std::vector< std::pair<unsigned int, double> > inRadiusIds;
	inRadiusIds.push_back( std::make_pair(id,0) );
	for ( unsigned int i = 0 ; i < radNeighbors.size() ; ++i )
	{
	    typename std::map<int, MeasurementVectorType>::iterator iter;
		MeasurementVectorType mvt = tree->GetMeasurementVector(radNeighbors[i]);
		for(iter = idToCentroidMap.begin(); iter != idToCentroidMap.end(); ++iter)
		{
			if(((*iter).second == mvt) && ((*iter).first!=id))
			{
				// to get the neighbors of all classes				
				if(Class_dest == 0)
				{
					double dist = distanceMetric->Evaluate( tree->GetMeasurementVector( radNeighbors[i] ));
					std::pair<unsigned int, double> Neighbor = std::make_pair((*iter).first, dist);
					inRadiusIds.push_back(Neighbor);
					break;
				}

				// to get the neighbors of a particular classes
				else
				{
					for(int col=((int)featureTable->GetNumberOfColumns())-1; col>=0; --col)
					{	
						std::string current_column = featureTable->GetColumnName(col);
						if(current_column.find("prediction") != std::string::npos )
						{
							for(int row = 0; row<(int)featureTable->GetNumberOfRows(); ++row)
							{
								if( ((*iter).first == featureTable->GetValue(row, 0).ToUnsignedInt())
									&& (featureTable->GetValue(row, col).ToUnsignedShort() == Class_dest) )
								{
									double dist = distanceMetric->Evaluate( tree->GetMeasurementVector( radNeighbors[i] ));
									std::pair<unsigned int, double> Neighbor = std::make_pair((*iter).first, dist);
									inRadiusIds.push_back(Neighbor);
									break;
								}
							}
							break;
						}
					}
					break;
				}
			}
		}
	}

	return inRadiusIds;
}




//***********************************************************************************************
//***********************************************************************************************
// Build GraphTable from Vectors
//***********************************************************************************************
//***********************************************************************************************
template <int num_dimensions>
vtkSmartPointer<vtkTable> kNearestObjects<num_dimensions>::vectorsToGraphTable(std::vector< std::vector< std::pair<unsigned int, double> > > NeighborIDs)
{
	graphtable = vtkSmartPointer<vtkTable>::New();
	graphtable->Initialize();
	vtkSmartPointer<vtkStringArray> column1 = vtkSmartPointer<vtkStringArray>::New();
    column1->SetName( "Source" );
    graphtable->AddColumn(column1);
    vtkSmartPointer<vtkStringArray> column2 = vtkSmartPointer<vtkStringArray>::New();
    column2->SetName( "Target" );
	graphtable->AddColumn(column2);
	vtkSmartPointer<vtkStringArray> column3 = vtkSmartPointer<vtkStringArray>::New();
	column3->SetName( "Distance" );
	graphtable->AddColumn(column3);

	allIds = 1;
	for(unsigned int j=0 ; j<NeighborIDs.size() ; j++)
	{
		this->vectorsToGraphTable(NeighborIDs[j]);
	}
	allIds = 0;

	return graphtable;
}

template <int num_dimensions>
vtkSmartPointer<vtkTable> kNearestObjects<num_dimensions>::vectorsToGraphTable(std::vector< std::pair<unsigned int, double> > NeighborIds)
{
	if(allIds == 0)
	{
		graphtable = vtkSmartPointer<vtkTable>::New();
		graphtable->Initialize();
		vtkSmartPointer<vtkStringArray> column1 = vtkSmartPointer<vtkStringArray>::New();
		column1->SetName( "Source" );
		graphtable->AddColumn(column1);
		vtkSmartPointer<vtkStringArray> column2 = vtkSmartPointer<vtkStringArray>::New();
		column2->SetName( "Target" );
		graphtable->AddColumn(column2);
		vtkSmartPointer<vtkStringArray> column3 = vtkSmartPointer<vtkStringArray>::New();
		column3->SetName( "Distance" );
		graphtable->AddColumn(column3);
	}
	
	for(unsigned int i=1 ; i<(int)NeighborIds.size() ; i++)
	{
		vtkVariant curSource = static_cast<vtkVariant>(NeighborIds.at(0).first);
		vtkVariant curTarget = static_cast<vtkVariant>(NeighborIds.at(i).first);
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
     	model_data1->InsertNextValue(curSource);
		model_data1->InsertNextValue(curTarget);
		model_data1->InsertNextValue(static_cast<vtkVariant>(NeighborIds.at(i).second));
		bool match = false;
		int j = 0;
		while ((j < this->graphtable->GetNumberOfRows())&&!match)
		{
			if ((graphtable->GetValue(j, 0) == curTarget)&&(graphtable->GetValue(j, 1) == curSource))
			{
				match = true;
			}
			else
			{
				j++;
			}
		}
		if (!match)
		{
			graphtable->InsertNextRow(model_data1);	
		}
	}
	return graphtable;
}




//***********************************************************************************************
//***********************************************************************************************
// Build Graph from Vectors
//***********************************************************************************************
//***********************************************************************************************
template <int num_dimensions>
typename kNearestObjects<num_dimensions>::NeighborGraph kNearestObjects<num_dimensions>::getNeighborGraph(std::vector< std::vector< std::pair<unsigned int, double> > > NeighborIDs)
{
	nodeName  = get(boost::vertex_name, NG);	
	allIds = 1;
	for(unsigned int j=0 ; j<NeighborIDs.size() ; j++)
	{
		this->getNeighborGraph(NeighborIDs[j]);
	}
	allIds = 0;

	return NG;
}

template <int num_dimensions>
typename kNearestObjects<num_dimensions>::NeighborGraph kNearestObjects<num_dimensions>::getNeighborGraph(std::vector< std::pair<unsigned int, double> > NeighborIds)
{
	if(allIds == 0)
	{
		NG.clear();
		nodeName  = get(boost::vertex_name, NG);
	}

	int src = GetNodeIndex(NeighborIds.at(0).first);
	if(src == -1)
	{
		node S;
		S = add_vertex(NG);
		nodeName[S] = convert2string(NeighborIds.at(0).first);
		src = num_vertices(NG)-1;
	}

	for(unsigned int i=1 ; i<(int)NeighborIds.size() ; i++)
	{
		int trg = GetNodeIndex(NeighborIds[i].first);
		if(trg ==-1)
		{
			node T;
			T = add_vertex(NG);
			nodeName[T] = convert2string(NeighborIds[i].first);
			trg = num_vertices(NG)-1;
		}
		bool bRet;
		Edge e;
		tie(e,bRet) = edge(src,trg,NG);
		if(!bRet)
		{
			add_edge(src,trg,NG);
		}
	}

	return NG;
}

template <int num_dimensions>
std::string kNearestObjects<num_dimensions>::convert2string(unsigned int id)
{
	std::stringstream out;
	out << id;
	std::string s = out.str();
	return s;
}

template <int num_dimensions>
int kNearestObjects<num_dimensions>::GetNodeIndex(unsigned int id)
{	
	std::string id1 = convert2string(id);
	int index;
	int flag = -1;
	for ( unsigned int i = 0 ; i < num_vertices(NG) ; ++i)
	{
		if(nodeName[i]==id1)
		{
			flag = 0;
			index = i;
			break;
		}
	}
	
	if (flag==0) 
	{
		return index;
	}
	else
	{
		return flag;
	}
}

//***********************************************************************************************
// Build Table from graph
//***********************************************************************************************
template <int num_dimensions>
vtkSmartPointer<vtkTable> kNearestObjects<num_dimensions>::graphToTable(NeighborGraph g)
{
	graphtable = vtkSmartPointer<vtkTable>::New();
	graphtable->Initialize();
	vtkSmartPointer<vtkStringArray> column1 = vtkSmartPointer<vtkStringArray>::New();
    column1->SetName( "Source" );
    graphtable->AddColumn(column1);
    vtkSmartPointer<vtkStringArray> column2 = vtkSmartPointer<vtkStringArray>::New();
    column2->SetName( "Target" );
	graphtable->AddColumn(column2);

	for (tie(vi,vi_end) = vertices(g) ; vi != vi_end ; ++vi)
    {
	    for (tie(ai,ai_end) = adjacent_vertices(*vi, g) ; ai != ai_end ; ++ai)
	    {
         //   check = 1;
	        //for(int i=0; i<(int)graphtable->GetNumberOfRows(); ++i)
	        //{
		       // if(((static_cast<vtkVariant>(nodeName[*vi])==graphtable->GetValue(i,0))&&(static_cast<vtkVariant>(nodeName[*ai])==graphtable->GetValue(i,1))) || 
		 	     //  ((static_cast<vtkVariant>(nodeName[*ai])==graphtable->GetValue(i,0))&&(static_cast<vtkVariant>(nodeName[*vi])==graphtable->GetValue(i,1))))
		       // {
		 	     //   check = 0;
		 	     //   break;
		       // }
 	       // }
	        //if(check == 1)
	        //{
	            vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
     	        model_data1->InsertNextValue(static_cast<vtkVariant>(nodeName[*vi]));
		        model_data1->InsertNextValue(static_cast<vtkVariant>(nodeName[*ai]));
		        graphtable->InsertNextRow(model_data1);		  
	        //}	   
        }
    }
    return graphtable;
}
	

#endif



