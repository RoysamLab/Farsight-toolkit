#include <iostream>
#include <string>

//#include "TraceXMLReader.h"
#include "View3D.h"
#include "Trace.h"
//#include "vtkPolyLine.h"

int main (int argc, char* argv[])	{
	View3d View; //initalizes the 3d viewer
	View.RenderWin();
  View.Initialize(argc, argv);
	vtkActor *act = View.LineAct();
	std::vector<TraceBit> vec = View.tobj->CollectTraceBits();
	printf("vec.size() = %d\n",vec.size());
//	View.AddPointsAsPoints(vec);
	View.AddBranchIllustrators();
	act->SetPickable(1);
	View.addAct(act);
	View.interact();
	printf("returned from View.interact\n");
	
	//tobject.WriteToSWCFile("../myswc.swc");
	
}
