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

#ifndef _TRACKSELECTION_H
#define _TRACKSELECTION_H
#include "helpers.h"
#include <set>
#include <QObject>


struct track_object{
	int tid;
	std::set<int> time_points;
};

struct lessthan{
	bool operator()(const track_object &a,const track_object &b) const
	{
		return a.tid<b.tid;
	}

};
/* Efficiently handles the selection of objects with time subsets */
class TrackSelection: public QObject{
	Q_OBJECT
public:
	TrackSelection(int n){m_selected.clear();max_time = n;}
	bool isSelected(int,int);
	bool add(int,int);
	bool remove(int,int);
	std::set<int> get(int);
	void setMaxTime(int);
	void clear(void);
signals:
	void changed(void);

private:
	int max_time;

	std::set<track_object,lessthan> m_selected;

};
#endif
