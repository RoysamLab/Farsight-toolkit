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

#include "TrackSelection.h"

//poorly implemented : May just need a sparse 2-D array
void TrackSelection::clear()
{
	m_selected.clear();
}
bool TrackSelection::isSelected(int num, int t=-1) // isSelected(.,-1) refers to "all selected" not "any selected"
{
	DEBUG3("isSelected: num = %d ",num);
	track_object tobj; tobj.tid = num;
	std::set<track_object,lessthan>::iterator iter = m_selected.find(tobj);
	if(iter==m_selected.end())
	{
		DEBUG2("returned false1\n");
		return false;
	}
	if(iter->time_points.size()==0)
	{
		DEBUG2("returned true1\n");
		return true;
	}
	if(t!=-1)
	{
		if(iter->time_points.find(t)==iter->time_points.end())
		{
			DEBUG2("returned false2\n");
			return false;
		}
		else
		{
			DEBUG2("returned true2\n");
			return true;
		}
	}
	else
	{
		DEBUG2("returned false3\n");
		return false;
	}
}

bool TrackSelection::add(int num, int t=-1)
{
	track_object tobj; tobj.tid = num;
	std::set<track_object,lessthan>::iterator iter = m_selected.find(tobj);
	if(iter==m_selected.end())
	{
		track_object tobj;
		tobj.tid = num;
		if(t!=-1)
		{
			tobj.time_points.insert(t);
		}
		m_selected.insert(tobj);
		emit changed();
		printf("1: returning true from add selection\n");
		return true;
	}
	if(t==-1)
	{
		tobj = *iter;
		m_selected.erase(iter);
		tobj.time_points.clear();
		m_selected.insert(tobj);
		printf("2: returning true from add selection\n");
	}
	else
	{
		tobj = *iter;
		m_selected.erase(iter);
		tobj.time_points.insert(t);
		if(tobj.time_points.size()==(max_time+1))
		{
			tobj.time_points.clear();
		}
		m_selected.insert(tobj);
	}
	emit changed();
	DEBUG3("Emitted changed()\n");
	return true;
}
bool TrackSelection::remove(int num, int t=-1)
{
	track_object tobj; tobj.tid = num;
	std::set<track_object,lessthan>::iterator iter = m_selected.find(tobj);
	if(iter==m_selected.end())
	{
		printf("Error: I should not be here in \"iter==m_selected.end()\" of TrackSelection::remove\n");
		return false;
	}
	if(t==-1)
	{
		m_selected.erase(iter);
	}
	else
	{
		tobj = *iter;
		m_selected.erase(iter);
		if(tobj.time_points.size()==0)
		{
			for(int counter=0; counter<=max_time; counter++)
			{
					tobj.time_points.insert(t);
			}
		}
		tobj.time_points.erase(t);
		if(tobj.time_points.size()!=0)
		{
			m_selected.insert(tobj);
		}
	}
	emit changed();
	DEBUG3("Emitted changed()\n");
	return true;
}

std::set<int> TrackSelection::get(int t=-1)
{
	std::set<track_object,lessthan>::iterator iter = m_selected.begin();
	std::set<int> output;
	while(iter!=m_selected.end())
	{
		if(t==-1)
		{
			output.insert(iter->tid);
			++iter;
			continue;
		}
		
		if(iter->time_points.size()==0)
		{
			output.insert(iter->tid);
		}
		else
		{
			if(iter->time_points.find(t)!=iter->time_points.end())
			{
				output.insert(iter->tid);
			}
		}
		++iter;
		continue;
	}
	return output;
}

void TrackSelection::setMaxTime(int tmax)
{
	max_time = tmax;
}