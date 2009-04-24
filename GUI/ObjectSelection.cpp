#include "ObjectSelection.h"


//Constructor
ObjectSelection::ObjectSelection(int time_max)
{
	m_selected.clear();
	t_min = 0;
	t_max = time_max;
}

void ObjectSelection::clear(void)
{
	m_selected.clear();
}

//If there are no t, or want to know if all time points do not set t
bool ObjectSelection::isSelected(int id, int t)
{
	if( t!=-1 && (t<t_min || t>t_max) )
		return false;

	//I create an instance of the object I'm looking for
	ObjectNode obj;
	obj.id = id;
	
	SelectedSetType::iterator iter = m_selected.find(obj);	//Try to find it in the container

	if( iter == m_selected.end() )			//Object is not in set
		return false;

	if( iter->timepoints.size() == 0 )		//Object is in set, and all time points are selected (so the one I am looking for is too)
		return true;
	
	if( t == -1 )							//Object is in set, not all time points are selected, but request was for all time points)
		return false;

	//Now I know I'm looking for a specific time point in one of the objects that does have some selections.
	//Lets see if there is an exact match:
	if( iter->timepoints.find(t) == iter->timepoints.end() )	//Time point was not found, return no match
		return false;
	else
		return true;											//Time point was found, return true

	//If I got here I can return false, I never found the desired object (should never get to this line
	return false;
}

//Leave t=-1 to select all time-points
bool ObjectSelection::add(int id, int t)
{
	if( t!=-1 && (t<t_min || t>t_max) )
		return false;

	//I create an instance of the object I'm hoping to add
	ObjectNode obj;
	obj.id = id;

	SelectedSetType::iterator iter = m_selected.find(obj);	//Try to find it in the container

	if( iter == m_selected.end() )	//Object is not in set, so add it!
	{
		if(t != -1)
			obj.timepoints.insert(t);

		m_selected.insert(obj);
		emit changed(id, t);
		return true;
	}
	
	//Found the object, now attempt to change its selections
	if( iter->timepoints.size() == 0 )		//All time points are already selected (so the one I want to set is too)
		return false;

	if( t == -1)							//Not all t's are selected, but I want them to be
	{
		iter->timepoints.clear();
		emit changed(id,t);
		return true;
	}

	//Not all t's are selected, and I only want to select one t, try to insert:
	iter->timepoints.insert(t);
	emit changed(id,t);
	return true;
}

//Leave t=-1 to remove all time-points
bool ObjectSelection::remove(int id, int t)
{
	if( t!=-1 && (t<t_min || t>t_max) )
		return false;

	//I create an instance of the object I'm hoping to add
	ObjectNode obj;
	obj.id = id;

	SelectedSetType::iterator iter = m_selected.find(obj);	//Try to find it in the container

	if( iter == m_selected.end() )	//Object is not in set, so cannot remove it!
		return false;

	if( t == -1)					//I want to remove the whole object
	{
		m_selected.erase(iter);
		emit changed(id,t);
		return true;
	}

	//Found the object, now remove the specific t requested
	if( iter->timepoints.size() == 0 )		//All time points are already selected.
	{
		//So I need to add all time points to the list (except the one I am removing)
		for( int i=t_min; i<=t_max; ++i )
		{
			if( i == t) continue;
			iter->timepoints.insert(i);
		}
		emit changed(id,t);
		return true;
	}

	//Just erase the t desired
	if( iter->timepoints.erase(t) )
	{
		emit changed(id,t);
		return true;
	}
	else
	{
		return false;	//t was not there to be erased
	}

	//SHOULD NEVER GET HERE
	return false;
}

//Returns all of the time points selected for this object
std::vector<int> ObjectSelection::getTimePoints(int object_id)	
{
	//I create an instance of the object I'm hoping to add
	ObjectNode obj;
	obj.id = object_id;

	SelectedSetType::iterator iter = m_selected.find(obj);	//Try to find it in the container

	//iterated through the whole array and did not find desired object:
	if( iter == m_selected.end() )
	{
		return std::vector<int>(0);
	}
	
	if ( iter->timepoints.size() == 0 )
	{
		//populate set with all time points and return
		std::vector<int> allT;
		for( int i=t_min; i<=t_max; ++i )
		{
			allT.push_back(i);
		}
		return allT;
	}

	//Copy the list to the vector:
	std::vector<int> rVect;
	std::set<int>::iterator t_iter;
	for(t_iter = iter->timepoints.begin(); t_iter != iter->timepoints.end(); ++t_iter )
	{
		rVect.push_back( *t_iter );
	}
	return rVect;
}

//Returns all of the objects for which all time points are selected
std::vector<int> ObjectSelection::getObjectsFull(void)		
{
	std::vector<int> objects;
	SelectedSetType::iterator iter;
	for(iter = m_selected.begin(); iter != m_selected.end(); ++iter)
	{
		if( iter->timepoints.size() == 0 )
			objects.push_back( iter->id );
	}
	return objects;
}

//Returns all of the objects for which only some of the time points are selected
std::vector<int> ObjectSelection::getObjectsPart(void)		
{
	std::vector<int> objects;
	SelectedSetType::iterator iter;
	for(iter = m_selected.begin(); iter != m_selected.end(); ++iter)
	{
		if( iter->timepoints.size() != 0 )
			objects.push_back( iter->id );
	}
	return objects;
}

//Leave t==-1 to return all objects that are selected for ANY t:
std::vector<int> ObjectSelection::getObjects(int t)
{
	std::vector<int> objects;
	SelectedSetType::iterator iter;
	for(iter=m_selected.begin(); iter != m_selected.end(); ++iter )
	{
		if( t == -1 )
		{
			objects.push_back(iter->id);
			continue;
		}

		//Need to look for specific t
		if(iter->timepoints.size() == 0)		//all t selected for this object
		{
			objects.push_back(iter->id);
		}
		else
		{
			if( iter->timepoints.find(t) != iter->timepoints.end() )	//object was found
			{
				objects.push_back(iter->id);
			}
		}
	}
	return objects;
}
