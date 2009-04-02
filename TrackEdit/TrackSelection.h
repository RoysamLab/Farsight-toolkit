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
