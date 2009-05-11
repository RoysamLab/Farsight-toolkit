#include "itkVector.h"
class TraceNode {
public:
	long ID;
	long TraceID;
	double radius;
	itk::Vector<double,3> loc;
	std::vector<long> nbrID;
	
	TraceNode()	{
		ID = -1;
		TraceID = -1;
		radius = 0.0;
		loc.Fill(0);
	}
};
