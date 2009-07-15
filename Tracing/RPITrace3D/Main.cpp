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

/////////////////////////////////////////////////////////////////////////////////
// FILE: main.cpp
//
// The main program of the tracking algorithm. 
// This project performs 3D tracking by following the strongest boundary
// signal not only over all valid directions and shift values (as before) but also 
// over all valid slices.
//
//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable STLport code

#include "Main.h"    // global variables for main only
//By Yousef: 10-25-2007
//The header files are needed for xml parsing
#include "tinyxml/tinyxml.h"
//#include <libxml/parser.h>
//#include <libxml/tree.h>
//#include <libxml/xmlreader.h>
#define MY_ENCODING "ISO-8859-1"
//#pragma comment( lib, "libxml2.lib" ) 
/////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// Function: main(int argc, char *argv[])
//
// The main function for the 3D tracing program
// 
// The Command line MUST follow the following format (order is not important)
// 
// -g <ConfigFileName> -I InputImageName
//
// where, InputImageName must be in BioRad PIC format.
//
// If the 3D image is represented as a sequence of 2D images each in PGM format, 
// then the user must specify the file name of the first image and the number of
// slices. For example, the following command line:
//
// -g <ConfigFileName> -I Input2DImage001.pgm 30
//
// will cause the program to trace a 3D image consisting of the slices
// InputImage001.pgm, .... , InputImage030.pgm
//
//

int main(int argc, char* argv[])
{
  string SummaryFileName;
  string fName;
  string fName2;
  int i;
  Timer OverallTime;

  cout << "RPITrace3D version-1.1  Aug 8th 2005" << endl << endl;

  if(argc < 2)
    {
    cerr << "Usage: RPITrace3D <config file>" << endl;
    exit(0);
    }
  //gConfig.ProcessCommandLine(argc,argv);
  gConfig.XmlReadSettings(argv[1]);

  ///////////////////////////////////////////////////////////////////////////////
  // Process the command line, read configuration file, and set config parameters
  ProcessCommandLine (argc, argv);

  // create an ofstream pointer of a log file
  std::string OutputPath = gConfig.GetOutputPath();
  std::string InputPath = gConfig.GetInputPath();
  std::string ImageName = gConfig.GetImageName();
  std::string ImageType = gConfig.GetImageType();
  fName = OutputPath + ImageName + "Log.txt";
  gLogFile.open (fName.c_str());  

  SummaryFileName = OutputPath + ImageName + "Summary.txt";

  /////////////////
  // 1. Initialization:
  //
  // Create the shift vectors, left and right templates All such objects
  // and variables are defined globally
  
  ///////////////////////
  // 3. The main image tracking cycle
  //
  // 3.1 Create a 3D image object from all image files
  // 3.2 Find seed points
  // 3.3 track the 3D image
  // 3.4 display the resulting image

  // Read all the images to be tracked into a single 3D image.
  std::string ImageFile = InputPath + ImageName + "." + ImageType;
  The3DImage = new C3DImage (ImageFile.c_str(), ImageType);   

  cout << "Initialization:" << endl;
  Initialize ();  
  int iSlices = The3DImage->m_iSlices;
  int iRows = The3DImage->m_iRows;
  int iCols = The3DImage->m_iCols;
  
  gpuchThe3DImageData = &(The3DImage->data[0][0][0]);
  glThe3DImageLength = iRows * iCols * iSlices;
  The3DImage->RemovePixels (5); 
  if( gConfig.GetWriteOutputFiles() )
    The3DImage->Histogram (OutputPath + ImageName + "_Histogram.txt");

  EstimateMedianAndStdDev (The3DImage->m_aiHistogram, gi3DMedian, gf3DStdDev);

  // Fill the extra slices added before the first slice and after the
  // the last slice
  The3DImage->FillPaddingSlices (gi3DMedian + 5 + 1);

  gLogFile << "3DMedian: " << gi3DMedian << " +/- " << gf3DStdDev << endl;

  /////////////////////////////////
  // Generate Projection Images

  cout << "\tGenerating Projection Images ... " << flush;

  ProjectionImage = new CImage(iRows, iCols);
  The3DImage->MaxProjectXY(*ProjectionImage);
  TrackImageXY = new CImage(*ProjectionImage);
  TrackImageXZ = new CImage(iSlices, iCols);
  TrackImageYZ = new CImage(iSlices, iRows);
  The3DImage->MaxProjectYZ(*TrackImageYZ);
  The3DImage->MaxProjectXZ(*TrackImageXZ);
  CanvasXY = new CImage(*TrackImageXY);
  CanvasXZ = new CImage(*TrackImageXZ);
  CanvasYZ = new CImage(*TrackImageYZ);
  //Draw_Projections();

  cout << "Done" << endl;

  ///////////////////////
  // 3.2 Find Seed Points
  // 
  // Such points will be used to start the tracking process. A seed point is also
  // assigned a direction. 

  // compute image statisitics
  cout << "\tComputing Image Statistics ... " << flush; 
  ProjectionImage->ComputeStatistics(giMedian, gfStdDev);

  gLogFile << "giMedian: " << giMedian << " +/- " << gfStdDev << endl;

  giMinSeedPixelValue = giMedian;// + gfStdDev;

  memset(gaiForegroundHistogram, 0, sizeof(int) * 256);
  memset(gaiBackgroundHistogram, 0, sizeof(int) * 256);

  cout << "Done" << endl;

  // Find the seeds for all the image slices. Notice that finding the seeds for
  // each slice, tracking it, and then finding the seeds for the next slice is 
  // better because it avoids all those already tracked vessels. However, we must
  // follow this approach so that we can set the threshold for the soma detection.
  cout << "\tSearching for Seed Points ... " << flush;

  int NumberOfUnverifiedSeeds = FindSeedPoints2(*TrackImageXY);
  cout << NumberOfUnverifiedSeeds << " seed candidates found" << endl;

  // an array of pointers to the seed points sorted according to their
  // gray level value
  gapSortedArrayOfSeedPoints = new CPoint * [iRows * iCols];
  cout << "\tGenerating statistics from seed candidates ... " << flush;
  VerifyAndSortAllSeedPoints();

  if ( cfg_mode_debug && gConfig.GetWriteOutputFiles() )
  {
    //Draw_PointsOnProjections(UnverifiedSeeds, "UnverifiedSeeds");
    //Draw_SeedCandidates(seed_candidates);
    //Draw_SeedPointsOnProjections();
    //Draw_PointsOnProjections(VerifiedSeedsCenter, "VerifiedSeeds");
    //Draw_SeedPointDirections();
  }

  ///////////////
  // Estimate Image Statistics
  // The foreground and bakground arrays were filled during the initial
  // seed verification process.
  EstimateMedianAndStdDev(gaiForegroundHistogram,
    giForegroundMedian,
    gfForegroundStdDev);
  EstimateMedianAndStdDev(gaiBackgroundHistogram,
    giBackgroundMedian,
    gfBackgroundStdDev);

  gLogFile << "Foreground: " << giForegroundMedian << " +/- "
    << gfForegroundStdDev << endl;
  gLogFile << "Background: " << giBackgroundMedian << " +/- "
    << gfBackgroundStdDev << endl;
  cout << "Done" << endl; 

  // Detect Soma ?  
  ////

  if ( gConfig.GetDetectSoma() )
  {   
    //LocateSomas3();         
    //Yousef: 01-27-2006  
    //try this
    SomaLabelsImage = new C3DImage(iSlices,iRows,iCols);
    LocateSomas3();
    //or this
    //By Yousef/////////////////////////////////////////////////////////////////
    //std::string SomaImageFile = "100upoint5%hippo25x1unmixed_for_render_microglia.pic";
    //Soma3DImage = new C3DImage (SomaImageFile.c_str(), ImageType);
    //SomaLabelsImage = new C3DImage(iSlices,iRows,iCols);
      //LocateSomas3_v4();
    ///////////////////////////////////////////////////////////////////////////
  }

  // each template column produces a response of 3*contrast.
  // if we assume that the background and the forground are separated
  // by 3*gfStdDev then we arrive at the following threshold for stopping
  // criteria.
  //  giSmallestResponseAccepted = (3.0*gfStdDev+0.5)*3*2*giUsedTemplateLength;
  //int iMinimumTemplateLength = gConfig.GetMinimumTemplateLength();
  //giSmallestResponseAccepted = 3 * 2 * iMinimumTemplateLength;

  giContrastThreshold = (int) (3 * gfStdDev + 0.5);

  if (giContrastThreshold < 2)
    giContrastThreshold = 2;

  ///////////////////////////////////////////////////////////////////////////
  // 3.3 Track the image 
  //
  // starting from a valid seed point, track the corresponding vessel and
  // and all vessels reachable from it. 
  // Each tracked vessel segment is represented as a list of cross_sections
  // a cross section is characterized by two boundary points, a center 
  // point and their associated directions.
  cout << "\nTracing the 3D image:" << endl;

  Timer trace_time;
  TraceA3DImage(*The3DImage); 

  //By Yousef
  //Detect Branch points
  if ( gConfig.GetDetectBranches() )
  {
    BranchPoints = new CBranches();
    BranchPoints->Search();
  }

  if(gTheVessels.m_iNumOfElements > 1 && !gConfig.GetDisableVesselMerging())
  {
    gTheVessels.ExtendVessels();  
    gTheVessels.BreakOnIntersectionPoints();

    //// if an intersection point involves two segments connected from their ends
    //// merge such segments
    gTheVessels.MergeVessels();
    gTheVessels.UpdateIntersectionPoints();
  }
  //gTheVessels.DeleteShortNetwork(giParam_Trace_MinSegmentLength);

  //gTheVessels.DeleteNarrowVessels(4);

  gLogFile << "Contrast Used: " << gfContrast << endl;
  gLogFile << "Tracing time: " << trace_time.elapsedSec() << " seconds"
    << endl;
  cout << "\tTracing time: " << trace_time.elapsedSec() << " mseconds"
    << endl;
  
  //if (cfg_mode_debug)
  //{
  //  Draw_BorderlineOnProjections("horizontal");
  //  Draw_BorderlineOnProjections("vertical");
  //}
  
  if( gConfig.GetWriteOutputFiles() ) {
    Draw_CenterlineOnProjections();
    if ( gConfig.GetDetectBranches() )
      Draw_BranchPoints();
  }

  if(gConfig.GetQA())
  {
    C3DImage * VesselnessImage;
    std::string vesselness_filename = gConfig.GetVesselnessImage();
    std::string vesselness_full_name = InputPath + vesselness_filename; 
    VesselnessImage = new C3DImage(vesselness_full_name,"pic");
    std::cout << "q: " << compute_Q(VesselnessImage,0.5) << std::endl;
  }

  //Draw_Residual3DImage();

  //if (cfg_mode_debug)
  //{
  //  Draw_Centerline3D();
  //  Draw_3DVessels();
  //}


  //Added by Yousef on 3-25-2009: Assign center points IDs and parent IDs. 
  //However, set the parent ID at the first point of each segment to -1 regardless of branching points
  //Also, if there is a branching point at the end of the vessel, reverse the IDs of its points
  
  //First, create a temporary image to hold the trace points
  int*** tmp_image = new int**[The3DImage->m_iSlices];
	for(i=0;i<The3DImage->m_iSlices;i++)
		tmp_image[i] = new int*[The3DImage->m_iRows];
	for(i=0;i<The3DImage->m_iSlices;i++)
	{
		for(int j=0;j<The3DImage->m_iRows;j++)
		{
			tmp_image[i][j] = new int[The3DImage->m_iCols];
			for(int k=0; k<The3DImage->m_iCols; k++)
				tmp_image[i][j][k] = 0;
		}
	}

  int Pnt_IDs = 0;
  for (i = 0; i < gTheVessels.m_iNumOfElements; i++)
  {
	  if (gTheVessels.m_apData[i])
      {     
		  if(gTheVessels.m_apData[i]->GetParentID()>0 && gTheVessels.m_apData[i]->GetParentLocation() !=0)
		  {
			  Pnt_IDs++;
			  CLNode<CPoint>* temp1 = gTheVessels.m_apData[i]->m_Center.tail; 
			  temp1->data->SetPointID(Pnt_IDs);			  
			  temp1->data->SetParentID(-1);
			  tmp_image[temp1->data->m_iZ][temp1->data->m_iY][temp1->data->m_iX] = Pnt_IDs;
			  temp1 = temp1->before;
			  while(temp1)
			  {
				  Pnt_IDs++;
				  temp1->data->SetPointID(Pnt_IDs);
				  temp1->data->SetParentID(Pnt_IDs-1);
				  tmp_image[temp1->data->m_iZ][temp1->data->m_iY][temp1->data->m_iX] = Pnt_IDs;
				  temp1 = temp1->before;
			  }
		  }
		  else
		  {
			  Pnt_IDs++;
			  CLNode<CPoint>* temp1 = gTheVessels.m_apData[i]->m_Center.head; 
			  temp1->data->SetPointID(Pnt_IDs);
			  tmp_image[temp1->data->m_iZ][temp1->data->m_iY][temp1->data->m_iX] = Pnt_IDs;
			  temp1->data->SetParentID(-1);
			  temp1 = temp1->after;
			  while(temp1)
			  {
				  Pnt_IDs++;
				  temp1->data->SetPointID(Pnt_IDs);
				  temp1->data->SetParentID(Pnt_IDs-1);
				  tmp_image[temp1->data->m_iZ][temp1->data->m_iY][temp1->data->m_iX] = Pnt_IDs;
				  temp1 = temp1->after;
			  }
		  }		  
	  }
  
  }
  //Added by Yousef on 3-25-2009: now, take care of the parent IDs at branch points
  for (i = 0; i < gTheVessels.m_iNumOfElements; i++)
  {
	  if (gTheVessels.m_apData[i] && gTheVessels.m_apData[i]->GetParentID()>0)
	  {
		  int PRID = tmp_image[gTheVessels.m_apData[i]->ParentBranchPoint->m_iZ][gTheVessels.m_apData[i]->ParentBranchPoint->m_iY][gTheVessels.m_apData[i]->ParentBranchPoint->m_iX];
			  
		  if(gTheVessels.m_apData[i]->GetParentLocation() ==0)
		  {
			  CLNode<CPoint>* temp1 = gTheVessels.m_apData[i]->m_Center.head; 			  
			  temp1->data->SetParentID(PRID);
		  }
		  else 
		  {
			  CLNode<CPoint>* temp1 = gTheVessels.m_apData[i]->m_Center.tail; 			  
			  temp1->data->SetParentID(PRID);
		  }		  
	  }
  }
  for(i=0;i<The3DImage->m_iSlices;i++)
	{
		for(int j=0;j<The3DImage->m_iRows;j++)
		{
			delete tmp_image[i][j];			
		}
		delete tmp_image[i];			
	}
    delete tmp_image;			


  if ( gConfig.GetWriteOutputFiles() ) 
  {
    /*fName = OutputPath + ImageName + "TracedPoints.txt";  
    ofstream points(fName.c_str());

    //spit out traced points
    for (i = 0; i < gTheVessels.m_iNumOfElements; i++)
    {
      if (gTheVessels.m_apData[i])
      {
        CLNode<CPoint>* temp = gTheVessels.m_apData[i]->m_Center.head;
        while (temp)
        {
          points << temp->data->m_iX << " " << temp->data->m_iY << " "
            << temp->data->m_iZ << " " 
            << temp->data->m_fHWidth << " " << temp->data->m_fVWidth << " "
            << (int)temp->data->m_iHDir << " " << (int)temp->data->m_iVDir << " "
            << temp->data->m_iValue << endl;
          temp = temp->after;
        }   
        points << -1 << " " << -1 << " "
          << -1 << endl;
      }
    }*/
    //By Yousef 10-25-2007: Print the results into an XML document
    /*
     * old libXML stuff:
    
    xmlDocPtr doc = NULL;       // document pointer 
    xmlNodePtr root_node = NULL, node = NULL, node2 = NULL; // node pointers 
    xmlDtdPtr dtd = NULL;       // DTD pointer

    LIBXML_TEST_VERSION;
       
    // Creates a new document, a node and set it as a root node
    */

    //new tinyxml stuff:

    /*TiXmlDocument doc;
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "", "");
    doc.LinkEndChild(decl);

    TiXmlElement *tracingOutput = new TiXmlElement("Trace");
    doc.LinkEndChild(tracingOutput);
    tracingOutput->SetAttribute("program", "RPITrace3D");
    tracingOutput->SetAttribute("version", "1.0");
    tracingOutput->SetAttribute("Number_of_Neurites", gTheVessels.m_iNumOfElements);

    TiXmlElement *traceLine;
    TiXmlElement *traceBit;
   
    for (i = 0; i < gTheVessels.m_iNumOfElements; i++)
      {		 
      if (gTheVessels.m_apData[i])
        {
        traceLine = new TiXmlElement("TraceLine");
        tracingOutput->LinkEndChild(traceLine);
        traceLine->SetAttribute("ID", i+1);
        traceLine->SetAttribute("length", gTheVessels.m_apData[i]->GetLength());
		traceLine->SetAttribute("Parent", gTheVessels.m_apData[i]->GetParentID());
        CLNode<CPoint>* temp = gTheVessels.m_apData[i]->m_Center.head;
        int IID = 0;
        while (temp)
          {
          IID++;
          traceBit = new TiXmlElement("TraceBit");
          traceLine->LinkEndChild(traceBit);
          traceBit->SetAttribute("ID", IID);
          traceBit->SetAttribute("X", temp->data->m_iX);
          traceBit->SetAttribute("Y", temp->data->m_iY);
          traceBit->SetAttribute("Z", temp->data->m_iZ);
          traceBit->SetAttribute("H_Direction", temp->data->m_iHDir);
          traceBit->SetAttribute("V_Direction", temp->data->m_iVDir);
          temp = temp->after;
          }           
        }
      }
    string filename = OutputPath + ImageName + "TracedPoints.xml";
    doc.SaveFile(filename.c_str());*/


	//Added by Yousef on 3-25-2009: write an swc file format
	fName = OutputPath + ImageName + "TracedPoints.swc";  
    ofstream points(fName.c_str());

    //spit out traced points
    for (i = 0; i < gTheVessels.m_iNumOfElements; i++)
    {
		if (gTheVessels.m_apData[i])
        {
			if(gTheVessels.m_apData[i]->GetParentLocation() == 0)
			{
				CLNode<CPoint>* temp = gTheVessels.m_apData[i]->m_Center.head;
				while (temp)
				{
					points << temp->data->GetPointID() << " " << 1 << " "
						   << temp->data->m_iX << " " << temp->data->m_iY << " " << temp->data->m_iZ << " " << 1 << " "
						   << temp->data->GetParentID() << " " <<endl;						   
					temp = temp->after;
				}   
			}
			else
			{
				CLNode<CPoint>* temp = gTheVessels.m_apData[i]->m_Center.tail;
				while (temp)
				{
					points << temp->data->GetPointID() << " " << 1 << " "
						   << temp->data->m_iX << " " << temp->data->m_iY << " " << temp->data->m_iZ << " " << 1 << " "
						   << temp->data->GetParentID() << " " <<endl;	
					temp = temp->after;
				}   
			}
		}
	}

    if ( gConfig.GetDetectSoma() )
    {
      if (!cfg_output_soma_draw_trees)
        gTheSomas->ConstructTrees();
      gTheSomas->Print(const_cast<char*> (SummaryFileName.c_str()));
    }
    else
    //  gTheVessels.Print(const_cast<char*> (SummaryFileName.c_str()));
    if(gConfig.GetDetectBranches())
    {
      fName2 = OutputPath + ImageName + "BranchPoints.txt"; 
      ofstream Br_points(fName2.c_str());
      
      for (i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
      {
        CPoint* pPoint =& gIntersectionPoints.m_apData[i]->m_Point;
        Br_points << pPoint->m_iX << " " << pPoint->m_iY << " " << pPoint->m_iZ <<endl;		
      }

    }

    cout << "Traces Summary: " << SummaryFileName << endl;
  }

  //try this
  //Draw_3DVessels();
  //Draw_Binary3DImage();
  ///
  
  FreeResources();
  

  cout << "Overall Execution Time: " << OverallTime.elapsedSec()
    << " seconds" << endl;
  gLogFile << "Overall Execution Time: " << OverallTime.elapsedSec()
    << " seconds" << endl;
  return 0;
}

