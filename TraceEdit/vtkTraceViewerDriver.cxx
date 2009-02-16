#include <iostream>
#include <string>


#include "View3D.h"
#include "Trace.h"
//#include "vtkPolyLine.h"



int main (int argc, char* argv[])	{
	//QApplication app(argc,argv);
	View3d View; //initalizes the 3d viewer
	TraceObject tobject;
	int num_loaded = 0;
	View.RenderWin();
	//tobject.ReadFromRPIXMLFile("../Yousef_TracedPoints.xml");
	for(int counter=1; counter<argc; counter++)// load as many files as possible. Provide offset for differentiating types
	{
		int len = strlen(argv[counter]);
		if(strcmp(argv[counter]+len-3,"swc")==0)
		{
			printf("I detected swc\n");
			tobject.ReadFromSWCFile(argv[counter]);
		}
		else if (strcmp(argv[counter]+len-3,"xml")==0)
		{
			printf("I detected xml\n");
			tobject.ReadFromRPIXMLFile(argv[counter]);
		}
		else if( strcmp(argv[counter]+len-3,"tks")==0)
		{
			printf("I detected tks\n");
			//tobject.ReadFromFeatureTracksFileForKymograph(argv[counter],num_loaded);
			tobject.ReadFromFeatureTracksFile(argv[counter],num_loaded);
		}
		else if( strcmp(argv[counter]+len-3,"tif")==0)
		{
			printf("I detected a tif file\n");
			//scanf("%*c");
			//View.readImg(argv[counter]);
			//View.AddContourThresholdSliders();
			View.rayCast(argv[counter]);
			View.AddVolumeSliders();
			//View.AddPlaybackWidget(argv[counter]);
			//scanf("%*c");
		}
		num_loaded++;
	}
	//tobject.ReadFromSWCFile("../n18fts.CNG.swc");
	//tobject.ReadFromSWCFile("../myswc.swc");
	//tobject.Print(std::cout);
	//TraceXMLReader::Pointer TraceReader = TraceXMLReader::New(); //creates xml reader
	//TraceReader->SetFileName("../traces.xml");
	//TraceReader->SetFileName(argv[1]);
	//TraceReader->Read();

		
	//vtkActor *act = View.LineAct(TraceReader->GetTraces());
	vtkActor *act = View.LineAct(tobject.GetVTKPolyData());
	std::vector<TraceBit> vec = tobject.CollectTraceBits();
	printf("vec.size() = %d\n",vec.size());
	View.AddPointsAsPoints(vec);
	View.tobj=&tobject;
	act->SetPickable(1);
	View.addAct(act);

	/*if(argc >= 3)
	{
		std::cout<<"\t source image set for contour filter \n";
		View.readImg(argv[2]);
	}
	if(argc >= 3)
	{
		std::cout<<"\t source image set for raycaster \n";
		View.rayCast(argv[2]);
	}*/
	View.interact();
	printf("returned from View.interact\n");
	//return app.exec();
	//tobject.WriteToSWCFile("../myswc.swc");
	
}
