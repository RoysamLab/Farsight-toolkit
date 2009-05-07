#ifndef _OBJECTSELECTION_H
#define _OBJECTSELECTION_H

#include <set>
#include <vector>
#include <QObject>

//****************************************************************************
// A class to keep track of selected objects/time points.
// If all of object is selected then timepoints will be empty
// use t=-1 to specify that all time points are desired (this is default)
// Assumes that timepoints in data are consecutive and begin with 0.
//**************************************************************************
class ObjectSelection: public QObject
{
	Q_OBJECT;

public:
	ObjectSelection(int time_max);
	
	void IncrementTMax(void);
	bool isSelected(int id, int t=-1);		
	bool add(int id, int t=-1);
	bool remove(int id, int t=-1);
	void clear(void);

	std::vector<int> getTimePoints(int object_id);	//Returns all of the time points selected for this object
	std::vector<int> getObjectsFull(void);			//Returns all of the objects IDs for which all time points are selected
	std::vector<int> getObjectsPart(void);			//Returns all of the objects IDs for which only some of the time points are selected
	std::vector<int> getObjects(int t=-1);			//Returns all objects that are selected for a given t value (default is objects selected at any t)

signals:
	void changed(int id, int t);
	
private:
	
	struct ObjectNode {
		int id;
		std::set<int> timepoints;
	};

	struct lessthan {
		bool operator()(const ObjectNode &a, const ObjectNode &b) const
		{
			return a.id < b.id;
		}
	};

	typedef std::set<ObjectNode,lessthan> SelectedSetType;

	SelectedSetType m_selected;
	int t_min;
	int t_max;

};


#endif