// SpineMain.cxx
#include "SpineRing.h"
//#include "itkMaximumProjectionImageFilter.h"
#include "ImageProc.h"

int main (int argc, char* argv[]) {
	const char *xmltrfilename = argv[2];
	
	//const char *xmltrfilename = "C:/FTK/spinering/detection/bin/debug/spine2bow.tif_SEs.txt.xml";
	ReaderType::Pointer reader = ReaderType::New();
//	WriterType::Pointer writer = WriterType::New();
   	reader->SetFileName( argv[1] );

	//reader->SetFileName( "C:/FTK/farsight/spinering/detection/bin/debug/spine2bow.tif");
	SpineImageType::Pointer image = reader->GetOutput();
	VERBOSE("Loading Image ...");
	image->Update();
	VERBOSE("Done loading image");
	VERBOSE("max int. z-proj ...");
	ImageDebugger imdbg(image, std::string(argv[1]));
	//imdbg.MIP = MIPGenerator<SpineImageType, SpineImage2DType>(image);
	//imdbg.MIP = MIPGenerator(image);
	TraceContainer TrCont(xmltrfilename, &imdbg);
	VERBOSE("finished reading xml file. Checking consistency..."); 
	TrCont.xmlconsistent = true;
	if (argc>3) {
		int chkconsist = atoi(argv[3]);
		if (chkconsist)
			TrCont.CheckConsistency(image);
	}
	//DetectorQ DQ; no Q's for now
	if (TrCont.xmlconsistent)    {
		GlobalDetector gd(image, &TrCont, &imdbg);
		VERBOSE("Detector running");
		gd.Run();
		
		VERBOSE("Validating Sp Cands");
		 
		gd.ValidateSpCands(); // & compute paths 
		imdbg.WriteSpCandidates2RGBImage(&gd);
		std::string pathfn = argv[1];
		pathfn += "RGBPaths.png";
		writeImageToFile(gd.GetRGBPathsImage(), pathfn);
 
//		gd.Clean();zzzzzzzzzzzzz
		//gd.SE(); later zzzzzzzzzzz
		//
//		gd.WriteXML();zzzzzzzzzzzz
	}

}