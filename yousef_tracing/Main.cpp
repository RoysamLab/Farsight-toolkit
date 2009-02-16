/////////////////////////////////////////////////////////////////////////////////
// FILE: main.cpp
//
// The main program of the tracking algorithm. 
// This project performs 3D tracking by following the strongest boundary
// signal not only over all valid directions and shift values (as before) but also 
// over all valid slices.
//
#pragma warning(disable:4786)  // disable STL-related warnings
#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
#pragma warning(disable:4702)  // unreachable STLport code

#include "Main.h"	   // global variables for main only

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
	int i;
	Timer overall_time;

	cout << "RPITrace3D version-1.1  June 20th 2005" << endl << endl;

	gConfig.ProcessCommandLine(argc,argv);

	///////////////////////////////////////////////////////////////////////////////
	// Process the command line, read configuration file, and set config parameters
	ProcessCommandLine (argc, argv);

	// create an ofstream pointer of a log file
	std::string output_path = gConfig.GetOutputPath();
	std::string input_path = gConfig.GetInputPath();
	std::string image_name = gConfig.GetImageName();
	std::string image_type = gConfig.GetImageType();
	fName = output_path + image_name + "Log.txt";
	gLogFile.open (fName.c_str());	

	SummaryFileName = output_path + image_name + "Summary.txt";

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
	std::string image_file = input_path + image_name + "." + image_type;
	The3DImage = new C3DImage (image_file.c_str(), image_type);

	cout << "Initialization:" << endl;
	Initialize ();	
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;

	gpuchThe3DImageData = &(The3DImage->data[0][0][0]);
	glThe3DImageLength = iRows * iCols * iSlices;
	The3DImage->RemovePixels (5);	
	if( gConfig.GetWriteOutputFiles() )
		The3DImage->Histogram (output_path + image_name + "_Histogram.txt");

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

	// rather than searching for seed points in every slice,
	// do so in the projection image first. Add all those verifiable seed points
	// to the array of seed points. 
	//std::vector<CPoint> seed_candidates;

	//if(cfg_seedsfromfile)
	//{
	//	candidates = ReadSeedPointsFromAFile(cfg_seedsfilename);
	//}
	//else
	//{
	//candidates = FindSeedPoints2(*TrackImageXY);
	FindSeedPoints2(*TrackImageXY);
	//}
	//int chunk_size = 15;
	//bool done = false;
	//int count = 1;
	//int z0, z1;
	//CImage * ImageForSeedDetection = new CImage(iRows, iCols);
	//while (!done)
	//{
	//	z0 = (chunk_size * count) - chunk_size;
	//	z1 = (chunk_size * count);
	//	
	//	The3DImage->MaxProjectXY(*ImageForSeedDetection,z0,z1);
	//	FindSeedPoints2(*ImageForSeedDetection, seed_candidates);
	//	std::cout << z0 << " to " << z1 << ": " << seed_candidates.size() << " found so far." << std::endl;
	//	if (z1 + chunk_size > iSlices)
	//		done = true;
	//	count++;
	//}
	//EstimateZCoord(seed_candidates);

	//std::vector<CPoint> XZ_seed_candidates;
	//count = 1;
	//int y0, y1;
	//CImage * XZ_ImageForSeedDetection = new CImage(iSlices,iCols);
	//done = false;
	//while (!done)
	//{
	//	y0 = (chunk_size * count) - chunk_size;
	//	y1 = (chunk_size * count);
	//	
	//	The3DImage->MaxProjectXZ(*XZ_ImageForSeedDetection,y0,y1);
	//	FindSeedPoints2(*XZ_ImageForSeedDetection, XZ_seed_candidates);
	//	std::cout << y0 << " to " << y1 << ": " << XZ_seed_candidates.size() << " found so far." << std::endl;
	//	if (y1 + chunk_size > iRows)
	//		done = true;
	//	count++;
	//}
	//for (i = 0; i < XZ_seed_candidates.size(); i++)
	//	XZ_seed_candidates[i].m_iZ = XZ_seed_candidates[i].m_iY;
	//EstimateYCoord(XZ_seed_candidates);

	//std::vector<CPoint> YZ_seed_candidates;
	//count = 1;
	//int x0, x1;
	//CImage * YZ_ImageForSeedDetection = new CImage(iSlices, iRows);
	//done = false;
	//while (!done)
	//{
	//	x0 = (chunk_size * count) - chunk_size;
	//	x1 = (chunk_size * count);
	//	
	//	The3DImage->MaxProjectYZ(*YZ_ImageForSeedDetection,x0,x1);
	//	FindSeedPoints2(*YZ_ImageForSeedDetection, YZ_seed_candidates);
	//	std::cout << x0 << " to " << x1 << ": " << YZ_seed_candidates.size() << " found so far." << std::endl;
	//	if (x1 + chunk_size > iCols)
	//		done = true;
	//	count++;
	//}
	//for (i = 0; i < YZ_seed_candidates.size(); i++)
	//{
	//	YZ_seed_candidates[i].m_iZ = YZ_seed_candidates[i].m_iY;
	//	YZ_seed_candidates[i].m_iY = YZ_seed_candidates[i].m_iX;
	//}
	//EstimateXCoord(YZ_seed_candidates);
	//
	//for (i = 0; i < XZ_seed_candidates.size(); i++)
	//	seed_candidates.push_back(XZ_seed_candidates[i]);
	//for (i = 0; i < YZ_seed_candidates.size(); i++)
	//	seed_candidates.push_back(YZ_seed_candidates[i]);

	cout << UnverifiedSeeds.size() << " seed candidates found" << endl;

	// an array of pointers to the seed points sorted according to thier
	// gray level value
	gapSortedArrayOfSeedPoints = new CPoint * [iRows * iCols];
	cout << "\tGenerating statistics from seed candidates ... " << flush;
	VerifyAndSortAllSeedPoints();
	//VerifyAndSortAllSeedPoints(seed_candidates);

	if (cfg_mode_debug && gConfig.GetWriteOutputFiles() )
	{
		Draw_PointsOnProjections(UnverifiedSeeds, "UnverifiedSeeds");
		//	cout << "DONE" << endl;
		//Draw_SeedCandidates(seed_candidates);
		Draw_SeedPointsOnProjections();
		Draw_PointsOnProjections(VerifiedSeedsCenter, "VerifiedSeeds");
		Draw_SeedPointDirections();
	}

	//============== seed clustering========================

	//CPoint seed1, seed2, seed3, seed4;
	//int maxX = -1; int minX = The3DImage->m_iCols;
	//int maxY = -1; int minY = The3DImage->m_iRows;
	//int number_of_seeds = VerifiedSeedsCenter.size();
	//int counter = 0;
	//int limit = (double)number_of_seeds;// / 3.0 ;
	//for (std::list<CPoint>::iterator s1 = VerifiedSeedsCenter.begin(); s1 != VerifiedSeedsCenter.end() ; s1++)
	//{
	//	if (s1->m_iX < minX) {
	//		minX = s1->m_iX;
	//		seed1 = *s1;
	//	}
	//	if (s1->m_iX > maxX) {
	//		maxX = s1->m_iX;
	//		seed2 = *s1;
	//	}
	//	if (s1->m_iY < minY) {
	//		minY = s1->m_iY;
	//		seed3 = *s1;
	//	}
	//	if (s1->m_iY > maxY) {
	//		maxY = s1->m_iY;
	//		seed4 = *s1;
	//	}
	//	counter++;
	//	if (counter > limit)
	//		break;
	//}

	//typedef itk::Vector<double, 3> MeasurementVectorType;
	//typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
	//SampleType::Pointer sample = SampleType::New();
	//MeasurementVectorType mv;
	//counter = 0;
	//for (std::list<CPoint>::iterator s1 = VerifiedSeedsCenter.begin(); s1 != VerifiedSeedsCenter.end() ; s1++)
	//{
	//	mv[0] = s1->m_iX;
	//	mv[1] = s1->m_iY;
	//	mv[2] = s1->m_iZ;
	//	sample->PushBack( mv );
	//	counter++;
	//	if (counter > limit)
	//		break;
	//}

	//typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
	//TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
	//treeGenerator->SetSample( sample );
	//treeGenerator->SetBucketSize( 16 );
	//treeGenerator->Update();

	//typedef TreeGeneratorType::KdTreeType TreeType;
	//typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
	//EstimatorType::Pointer estimator = EstimatorType::New();
	//EstimatorType::ParametersType initialMeans(12);
	//initialMeans[0] = seed1.m_iX;
	//initialMeans[1] = seed1.m_iY;
	//initialMeans[2] = seed1.m_iZ;
	//initialMeans[3] = seed2.m_iX;
	//initialMeans[4] = seed2.m_iY;
	//initialMeans[5] = seed2.m_iZ;
	//initialMeans[6] = seed3.m_iX;
	//initialMeans[7] = seed3.m_iY;
	//initialMeans[8] = seed3.m_iZ;
	//initialMeans[9] = seed4.m_iX;
	//initialMeans[10] = seed4.m_iY;
	//initialMeans[11] = seed4.m_iZ;

	//estimator->SetParameters( initialMeans );
	//estimator->SetKdTree( treeGenerator->GetOutput() );
	//estimator->SetMaximumIteration( 200 );
	//estimator->SetCentroidPositionChangesThreshold(0.0);
	//estimator->StartOptimization();

	//EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
	//for ( unsigned int i = 0 ; i < 12 ; ++i )
	//{
	//	std::cout << "cluster[" << i << "] " << std::endl;
	//	std::cout << " estimated mean : " << estimatedMeans[i] << std::endl;
	//}

	//std::ofstream out3("regionmeans.txt");
	//int i8 = 0;
	//out3 << estimatedMeans[i8++] << " "; out3 << estimatedMeans[i8++] << " "; out3 << estimatedMeans[i8++] << std::endl;
	//out3 << estimatedMeans[i8++] << " "; out3 << estimatedMeans[i8++] << " "; out3 << estimatedMeans[i8++] << std::endl;
	//out3 << estimatedMeans[i8++] << " "; out3 << estimatedMeans[i8++] << " "; out3 << estimatedMeans[i8++] << std::endl;
	//out3 << estimatedMeans[i8++] << " "; out3 << estimatedMeans[i8++] << " "; out3 << estimatedMeans[i8++] << std::endl;

	//typedef itk::Statistics::EuclideanDistance< MeasurementVectorType >
	//	MembershipFunctionType;
	//typedef itk::MinimumDecisionRule DecisionRuleType;
	//DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
	//typedef itk::Statistics::SampleClassifier< SampleType > ClassifierType;
	//ClassifierType::Pointer classifier = ClassifierType::New();
	//classifier->SetDecisionRule( (itk::DecisionRuleBase::Pointer) decisionRule);
	//classifier->SetSample( sample );
	//classifier->SetNumberOfClasses( 4 );
	//std::vector< unsigned int > classLabels;
	//classLabels.resize( 4 );
	//classLabels[0] = 1;
	//classLabels[1] = 2;
	//classLabels[2] = 3;
	//classLabels[3] = 4;
	//classifier->SetMembershipFunctionClassLabels( classLabels );

	//std::vector< MembershipFunctionType::Pointer > membershipFunctions;
	//MembershipFunctionType::OriginType origin;
	//int index = 0;
	//for ( unsigned int i = 0 ; i < 4 ; i++ )
	//{
	//	membershipFunctions.push_back( MembershipFunctionType::New() );
	//	for ( unsigned int j = 0 ; j < SampleType::MeasurementVectorSize ; j++ )
	//	{
	//		origin[j] = estimatedMeans[index++];
	//	}
	//	membershipFunctions[i]->SetOrigin( origin );
	//	classifier->AddMembershipFunction( membershipFunctions[i].GetPointer() );
	//}
	//classifier->Update();
	//ClassifierType::OutputType* membershipSample = classifier->GetOutput();
	//ClassifierType::OutputType::ConstIterator iter = membershipSample->Begin();
	//std::ofstream seed_clusters_out("seed_clusters.txt");
	//std::vector<int> count(4,0);

	//int X1_min = The3DImage->m_iCols; int X1_max = -1;
	//int Y1_min = The3DImage->m_iRows; int Y1_max = -1;
	//int Z1_min = The3DImage->m_iSlices; int Z1_max = -1;
	//int X2_min = The3DImage->m_iCols; int X2_max = -1;
	//int Y2_min = The3DImage->m_iRows; int Y2_max = -1;
	//int Z2_min = The3DImage->m_iSlices; int Z2_max = -1;
	//int X3_min = The3DImage->m_iCols; int X3_max = -1;
	//int Y3_min = The3DImage->m_iRows; int Y3_max = -1;
	//int Z3_min = The3DImage->m_iSlices; int Z3_max = -1;
	//int X4_min = The3DImage->m_iCols; int X4_max = -1;
	//int Y4_min = The3DImage->m_iRows; int Y4_max = -1;
	//int Z4_min = The3DImage->m_iSlices; int Z4_max = -1;
	//
	//typedef itk::Vector< float, 6 > RegionMeasurementVectorType;
	//typedef itk::Statistics::ListSample< RegionMeasurementVectorType > RegionSampleType;
	//RegionSampleType::Pointer R1_samples = RegionSampleType::New();
	//RegionSampleType::Pointer R2_samples = RegionSampleType::New();
	//RegionSampleType::Pointer R3_samples = RegionSampleType::New();
	//RegionSampleType::Pointer R4_samples = RegionSampleType::New();
	//RegionMeasurementVectorType rmv;
	//while ( iter != membershipSample->End() )
	//{
	//	seed_clusters_out << "measurement vector = " << iter.GetMeasurementVector()
	//		<< "class label = " << iter.GetClassLabel()
	//		<< std::endl;
	//	count[iter.GetClassLabel()-1]++;

	//	switch (iter.GetClassLabel()) 
	//	{
	//	case 1:
	//		X1_min = (iter.GetMeasurementVector()[0] < X1_min) ? iter.GetMeasurementVector()[0] : X1_min;
	//		X1_max = (iter.GetMeasurementVector()[0] > X1_max) ? iter.GetMeasurementVector()[0] : X1_max;
	//		Y1_min = (iter.GetMeasurementVector()[1] < Y1_min) ? iter.GetMeasurementVector()[1] : Y1_min;
	//		Y1_max = (iter.GetMeasurementVector()[1] > Y1_max) ? iter.GetMeasurementVector()[1] : Y1_max;
	//		Z1_min = (iter.GetMeasurementVector()[2] < Z1_min) ? iter.GetMeasurementVector()[2] : Z1_min;
	//		Z1_max = (iter.GetMeasurementVector()[2] > Z1_max) ? iter.GetMeasurementVector()[2] : Z1_max;
	//		for (std::list<CPoint>::iterator s1 = VerifiedSeedsCenter.begin(); s1 != VerifiedSeedsCenter.end() ; s1++)
	//			if( iter.GetMeasurementVector()[0] == s1->m_iX &&
	//				iter.GetMeasurementVector()[1] == s1->m_iY &&
	//				iter.GetMeasurementVector()[2] == s1->m_iZ )
	//			{ 
	//				rmv[0] = (s1->m_iHDir > 16) ? s1->m_iHDir-16 : s1->m_iHDir;
	//				rmv[1] = (s1->m_iVDir > 16) ? s1->m_iVDir-16 : s1->m_iVDir;
	//				rmv[2] = s1->m_fHWidth;
	//				rmv[3] = s1->m_fVWidth;
	//				rmv[4] = s1->m_iValue;
	//				rmv[5] = s1->m_iValue;
	//				R1_samples->PushBack( rmv );
	//			}
	//			break;
	//	case 2:
	//		X2_min = (iter.GetMeasurementVector()[0] < X2_min) ? iter.GetMeasurementVector()[0] : X2_min;
	//		X2_max = (iter.GetMeasurementVector()[0] > X2_max) ? iter.GetMeasurementVector()[0] : X2_max;
	//		Y2_min = (iter.GetMeasurementVector()[1] < Y2_min) ? iter.GetMeasurementVector()[1] : Y2_min;
	//		Y2_max = (iter.GetMeasurementVector()[1] > Y2_max) ? iter.GetMeasurementVector()[1] : Y2_max;
	//		Z2_min = (iter.GetMeasurementVector()[2] < Z2_min) ? iter.GetMeasurementVector()[2] : Z2_min;
	//		Z2_max = (iter.GetMeasurementVector()[2] > Z2_max) ? iter.GetMeasurementVector()[2] : Z2_max;
	//		for (std::list<CPoint>::iterator s1 = VerifiedSeedsCenter.begin(); s1 != VerifiedSeedsCenter.end() ; s1++)
	//			if( iter.GetMeasurementVector()[0] == s1->m_iX &&
	//				iter.GetMeasurementVector()[1] == s1->m_iY &&
	//				iter.GetMeasurementVector()[2] == s1->m_iZ )
	//			{ 
	//				rmv[0] = (s1->m_iHDir > 16) ? s1->m_iHDir-16 : s1->m_iHDir;
	//				rmv[1] = (s1->m_iVDir > 16) ? s1->m_iVDir-16 : s1->m_iVDir;
	//				rmv[2] = s1->m_fHWidth;
	//				rmv[3] = s1->m_fVWidth;
	//				rmv[4] = s1->m_iValue;
	//				rmv[5] = s1->m_iValue;
	//				R2_samples->PushBack( rmv );
	//			}
	//			break;
	//	case 3:
	//		X3_min = (iter.GetMeasurementVector()[0] < X3_min) ? iter.GetMeasurementVector()[0] : X3_min;
	//		X3_max = (iter.GetMeasurementVector()[0] > X3_max) ? iter.GetMeasurementVector()[0] : X3_max;
	//		Y3_min = (iter.GetMeasurementVector()[1] < Y3_min) ? iter.GetMeasurementVector()[1] : Y3_min;
	//		Y3_max = (iter.GetMeasurementVector()[1] > Y3_max) ? iter.GetMeasurementVector()[1] : Y3_max;
	//		Z3_min = (iter.GetMeasurementVector()[2] < Z3_min) ? iter.GetMeasurementVector()[2] : Z3_min;
	//		Z3_max = (iter.GetMeasurementVector()[2] > Z3_max) ? iter.GetMeasurementVector()[2] : Z3_max;
	//		for (std::list<CPoint>::iterator s1 = VerifiedSeedsCenter.begin(); s1 != VerifiedSeedsCenter.end() ; s1++)
	//			if( iter.GetMeasurementVector()[0] == s1->m_iX &&
	//				iter.GetMeasurementVector()[1] == s1->m_iY &&
	//				iter.GetMeasurementVector()[2] == s1->m_iZ )
	//			{ 
	//				rmv[0] = (s1->m_iHDir > 16) ? s1->m_iHDir-16 : s1->m_iHDir;
	//				rmv[1] = (s1->m_iVDir > 16) ? s1->m_iVDir-16 : s1->m_iVDir;
	//				rmv[2] = s1->m_fHWidth;
	//				rmv[3] = s1->m_fVWidth;
	//				rmv[4] = s1->m_iValue;
	//				rmv[5] = s1->m_iValue;
	//				R3_samples->PushBack( rmv );
	//			}
	//			break;
	//	case 4:
	//		X4_min = (iter.GetMeasurementVector()[0] < X4_min) ? iter.GetMeasurementVector()[0] : X4_min;
	//		X4_max = (iter.GetMeasurementVector()[0] > X4_max) ? iter.GetMeasurementVector()[0] : X4_max;
	//		Y4_min = (iter.GetMeasurementVector()[1] < Y4_min) ? iter.GetMeasurementVector()[1] : Y4_min;
	//		Y4_max = (iter.GetMeasurementVector()[1] > Y4_max) ? iter.GetMeasurementVector()[1] : Y4_max;
	//		Z4_min = (iter.GetMeasurementVector()[2] < Z4_min) ? iter.GetMeasurementVector()[2] : Z4_min;
	//		Z4_max = (iter.GetMeasurementVector()[2] > Z4_max) ? iter.GetMeasurementVector()[2] : Z4_max;
	//		for (std::list<CPoint>::iterator s1 = VerifiedSeedsCenter.begin(); s1 != VerifiedSeedsCenter.end() ; s1++)
	//			if( iter.GetMeasurementVector()[0] == s1->m_iX &&
	//				iter.GetMeasurementVector()[1] == s1->m_iY &&
	//				iter.GetMeasurementVector()[2] == s1->m_iZ )
	//			{ 
	//				rmv[0] = (s1->m_iHDir > 16) ? s1->m_iHDir-16 : s1->m_iHDir;
	//				rmv[1] = (s1->m_iVDir > 16) ? s1->m_iVDir-16 : s1->m_iVDir;
	//				rmv[2] = s1->m_fHWidth;
	//				rmv[3] = s1->m_fVWidth;
	//				rmv[4] = s1->m_iValue;
	//				rmv[5] = s1->m_iValue;
	//				R4_samples->PushBack( rmv );
	//			}
	//		break;			
	//	}
	//	++iter;
	//}

	//typedef itk::Statistics::MeanCalculator< RegionSampleType > MeanCalculatorType;
	//MeanCalculatorType::Pointer R1_mean = MeanCalculatorType::New();
	//MeanCalculatorType::Pointer R2_mean = MeanCalculatorType::New();
	//MeanCalculatorType::Pointer R3_mean = MeanCalculatorType::New();
	//MeanCalculatorType::Pointer R4_mean = MeanCalculatorType::New();
	//R1_mean->SetInputSample( R1_samples );
	//R2_mean->SetInputSample( R2_samples );
	//R3_mean->SetInputSample( R3_samples );
	//R4_mean->SetInputSample( R4_samples );
	//R1_mean->Update();
	//R2_mean->Update();
	//R3_mean->Update();
	//R4_mean->Update();
	//std::cout << "R1 mean: " << std::endl;
	//std::cout << *(R1_mean->GetOutput()) << std::endl;
	//std::cout << "R2 mean: " << std::endl;
	//std::cout << *(R2_mean->GetOutput()) << std::endl;
	//std::cout << "R3 mean: " << std::endl;
	//std::cout << *(R3_mean->GetOutput()) << std::endl;
	//std::cout << "R4 mean: " << std::endl;
	//std::cout << *(R4_mean->GetOutput()) << std::endl;

	//
	//typedef itk::Statistics::CovarianceCalculator< RegionSampleType > CovarianceCalculatorType;
	//CovarianceCalculatorType::Pointer R1_cov = CovarianceCalculatorType::New();
	//CovarianceCalculatorType::Pointer R2_cov = CovarianceCalculatorType::New();
	//CovarianceCalculatorType::Pointer R3_cov = CovarianceCalculatorType::New();
	//CovarianceCalculatorType::Pointer R4_cov = CovarianceCalculatorType::New();
	//R1_cov->SetInputSample( R1_samples );
	//R2_cov->SetInputSample( R2_samples );
	//R3_cov->SetInputSample( R3_samples );
	//R4_cov->SetInputSample( R4_samples );
	//R1_cov->SetMean( R1_mean->GetOutput() );
	//R2_cov->SetMean( R2_mean->GetOutput() );
	//R3_cov->SetMean( R3_mean->GetOutput() );
	//R4_cov->SetMean( R4_mean->GetOutput() );
	//R1_cov->Update();
	//R2_cov->Update();
	//R3_cov->Update();
	//R4_cov->Update();
	//
	//std::cout << "R1 cov: " << std::endl;
	//std::cout << *(R1_cov->GetOutput()) << std::endl;
	//std::cout << "R2 cov: " << std::endl;
	//std::cout << *(R2_cov->GetOutput()) << std::endl;
	//std::cout << "R3 cov: " << std::endl;
	//std::cout << *(R3_cov->GetOutput()) << std::endl;
	//std::cout << "R4 cov: " << std::endl;
	//std::cout << *(R4_cov->GetOutput()) << std::endl;

	//C3DImage * VImage;
	//std::string v_filename = gConfig.GetVesselnessImage();
	//std::string v_full_name = input_path + v_filename; 
	//VImage = new C3DImage(v_full_name,"pic");
	//int iPadding = gConfig.GetImagePadding();
	//// TODO: take account the boundary seeds' widths 

	//// create the volumes
	//int R1_cols = X1_max - X1_min;
	//int R1_rows = Y1_max - Y1_min;
	//int R1_slices = Z1_max - Z1_min;
	//C3DImage R1(R1_slices, R1_rows, R1_cols);
	//C3DImage R1V(R1_slices, R1_rows, R1_cols);
	//for (int s = Z1_min; s < Z1_max; s++){
	//	for (int r = Y1_min; r < Y1_max; r++) {
	//		for (int c = X1_min; c < X1_max; c++) {
	//			R1.data[ s - Z1_min ][r - Y1_min][c - X1_min] = The3DImage->data[s][r][c];
	//			R1V.data[ s - Z1_min ][r - Y1_min][c - X1_min] = VImage->data[s][r][c];
	//		}
	//	}
	//}
	//R1.Write("R1.pic");
	//R1V.Write("R1V.pic");

	//int R2_cols = X2_max - X2_min;
	//int R2_rows = Y2_max - Y2_min;
	//int R2_slices = Z2_max - Z2_min;
	//C3DImage R2(R2_slices, R2_rows, R2_cols);
	//C3DImage R2V(R2_slices, R2_rows, R2_cols);
	//for (int s = Z2_min; s < Z2_max; s++){
	//	for (int r = Y2_min; r < Y2_max; r++) {
	//		for (int c = X2_min; c < X2_max; c++) {
	//			R2.data[ s - Z2_min ][r - Y2_min][c - X2_min] = The3DImage->data[s][r][c];
	//			R2V.data[ s - Z2_min ][r - Y2_min][c - X2_min] = VImage->data[s][r][c];
	//		}
	//	}
	//}
	//R2.Write("R2.pic");
	//R2V.Write("R2V.pic");

	//int R3_cols = X3_max - X3_min;
	//int R3_rows = Y3_max - Y3_min;
	//int R3_slices = Z3_max - Z3_min;
	//C3DImage R3(R3_slices, R3_rows, R3_cols);
	//C3DImage R3V(R3_slices, R3_rows, R3_cols);
	//for (int s = Z3_min; s < Z3_max; s++){
	//	for (int r = Y3_min; r < Y3_max; r++) {
	//		for (int c = X3_min; c < X3_max; c++) {
	//			R3.data[ s - Z3_min ][r - Y3_min][c - X3_min] = The3DImage->data[s][r][c];
	//			R3V.data[ s - Z3_min ][r - Y3_min][c - X3_min] = VImage->data[s][r][c];
	//		}
	//	}
	//}
	//R3.Write("R3.pic");
	//R3V.Write("R3V.pic");

	//int R4_cols = X4_max - X4_min;
	//int R4_rows = Y4_max - Y4_min;
	//int R4_slices = Z4_max - Z4_min;
	//C3DImage R4(R4_slices, R4_rows, R4_cols);
	//C3DImage R4V(R4_slices, R4_rows, R4_cols);
	//for (int s = Z4_min; s < Z4_max; s++){
	//	for (int r = Y4_min; r < Y4_max; r++) {
	//		for (int c = X4_min; c < X4_max; c++) {
	//			R4.data[ s - Z4_min ][r - Y4_min][c - X4_min] = The3DImage->data[s][r][c];
	//			R4V.data[ s - Z4_min ][r - Y4_min][c - X4_min] = VImage->data[s][r][c];
	//		}
	//	}
	//}
	//R4.Write("R4.pic");
	//R4V.Write("R4V.pic");

	//CImage * RegionImage = new CImage(*ProjectionImage);
	//for( int x = X1_min; x <= X1_max; x++)
	//	RegionImage->data[Y1_min][x] = 0;
	//for( int x = X1_min; x <= X1_max; x++)
	//	RegionImage->data[Y1_max][x] = 0;
	//for( int y = Y1_min; y <= Y1_max; y++)
	//	RegionImage->data[y][X1_min] = 0;
	//for( int y = Y1_min; y <= Y1_max; y++)
	//	RegionImage->data[y][X1_max] = 0;
	//for( int x = X2_min; x <= X2_max; x++)
	//	RegionImage->data[Y2_min][x] = 1;
	//for( int x = X2_min; x <= X2_max; x++)
	//	RegionImage->data[Y2_max][x] = 1;
	//for( int y = Y2_min; y <= Y2_max; y++)
	//	RegionImage->data[y][X2_min] = 1;
	//for( int y = Y2_min; y <= Y2_max; y++)
	//	RegionImage->data[y][X2_max] = 1;
	//for( int x = X3_min; x <= X3_max; x++)
	//	RegionImage->data[Y3_min][x] = 2;
	//for( int x = X3_min; x <= X3_max; x++)
	//	RegionImage->data[Y3_max][x] = 2;
	//for( int y = Y3_min; y <= Y3_max; y++)
	//	RegionImage->data[y][X3_min] = 2;
	//for( int y = Y3_min; y <= Y3_max; y++)
	//	RegionImage->data[y][X3_max] = 2;
	//for( int x = X4_min; x <= X4_max; x++)
	//	RegionImage->data[Y4_min][x] = 3;
	//for( int x = X4_min; x <= X4_max; x++)
	//	RegionImage->data[Y4_max][x] = 3;
	//for( int y = Y4_min; y <= Y4_max; y++)
	//	RegionImage->data[y][X4_min] = 3;
	//for( int y = Y4_min; y <= Y4_max; y++)
	//	RegionImage->data[y][X4_max] = 3;
	//RegionImage->WriteTIFF("RegionImage.tif");
	//

	//for (int i = 0; i < 4; i++)
	//{
	//	std::cout << "Cluster " << i + 1 << " count: " << count[i] << std::endl;
	//}

	//======================================================

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
	//Yousef: 01-27-2006
	//SomaLabelsImage = new C3DImage(iSlices,iRows,iCols);
	////
	//if ( gConfig.GetDetectSoma() )
	//{
		//cout << "\tDetecting Somas ... ";
		//LocateSomas3();			
		//try this
		//LocateSomas3_v3();		
		//or this
		//By Yousef/////////////////////////////////////////////////////////////////
		//std::string SomaImageFile = "l100um2percent25x3unmixed_for_render_Iba1.pic";
		//std::string SomaImageFile = input_path + "100um2percent25x3unmixed_for_render_Iba1.pic";		
		//Soma3DImage = new C3DImage (SomaImageFile.c_str(), image_type);
		//SomaLabelsImage = new C3DImage(iSlices,iRows,iCols);
		//LocateSomas3_v4();
		///////////////////////////////////////////////////////////////////////////
	//}

	// each template column produces a response of 3*contrast.
	// if we assume that the background and the forground are separated
	// by 3*gfStdDev then we arrive at the following threshold for stopping
	// criteria.
	//	giSmallestResponseAccepted = (3.0*gfStdDev+0.5)*3*2*giUsedTemplateLength;
	int iMinimumTemplateLength = gConfig.GetMinimumTemplateLength();
	giSmallestResponseAccepted = 3 * 2 * iMinimumTemplateLength;

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
	//BranchPoints = new CBranches();
	//BranchPoints->Search();

	if(gTheVessels.m_iNumOfElements > 1 && !gConfig.GetDisableVesselMerging())
	{
		gTheVessels.ExtendVessels();	
		gTheVessels.BreakOnIntersectionPoints();

		//// if an intersection point involves two segments connected from their ends
		//// merge such segments
		gTheVessels.MergeVessels();
		gTheVessels.UpdateIntersectionPoints();
	}
	gTheVessels.DeleteShortNetwork(giParam_Trace_MinSegmentLength);

	//gTheVessels.DeleteNarrowVessels(4);

	gLogFile << "Contrast Used: " << gfContrast << endl;
	gLogFile << "Tracing time: " << trace_time.elapsedSec() << " seconds"
		<< endl;
	cout << "\tTracing time: " << trace_time.elapsedSec() << " mseconds"
		<< endl;

	//if (cfg_mode_debug)
	//{
	//	Draw_BorderlineOnProjections("horizontal");
	//	Draw_BorderlineOnProjections("vertical");
	//}
	if( gConfig.GetWriteOutputFiles() ) 
	{
		Draw_CenterlineOnProjections();
		//Draw_BranchPoints();
	}

	if(gConfig.GetQA())
	{
		C3DImage * VesselnessImage;
		std::string vesselness_filename = gConfig.GetVesselnessImage();
		std::string vesselness_full_name = input_path + vesselness_filename; 
		VesselnessImage = new C3DImage(vesselness_full_name,"pic");
		std::cout << "q: " << compute_Q(VesselnessImage,0.5) << std::endl;
	}

	//Draw_Residual3DImage();

	//if (cfg_mode_debug)
	//{
	//	Draw_Centerline3D();
	//	Draw_3DVessels();
	//}

	if ( gConfig.GetWriteOutputFiles() ) 
	{
		fName = output_path + image_name + "TracedPoints.txt";	
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
		}
		
		/*if ( gConfig.GetDetectSoma() )
		{
			if (!cfg_output_soma_draw_trees)
				gTheSomas->ConstructTrees();
			gTheSomas->Print(const_cast<char*> (SummaryFileName.c_str()));
		}
		else*/
			gTheVessels.Print(const_cast<char*> (SummaryFileName.c_str()));
		
		cout << "Traces Summary: " << SummaryFileName << endl;
		
		//for (i = 0; i < gTheVessels.m_iNumOfElements; i++)
		//{
		//	int length = 0;
		//	if (gTheVessels.m_apData[i])
		//	{
		//		CLNode<CPoint>* temp = gTheVessels.m_apData[i]->m_Center.head;
		//		while (temp)
		//		{
		//			length++;
		//			temp = temp->after;
		//		}
		//		points << length << std::endl;
		//		temp = gTheVessels.m_apData[i]->m_Center.head;
		//		while (temp)
		//		{
		//			points << temp->data->m_iX << " " << temp->data->m_iY << " "
		//				<< temp->data->m_iZ << " " << temp->data->m_fHWidth << " " << temp->data->m_fVWidth << endl;
		//			temp = temp->after;
		//		}
		//	}
		//}

		//if (cfg_mode_debug)
		//{
		//	Draw_PointsOnProjections(CenterPoints, HLPoints, HRPoints, "trace");
		//	Draw_PointsOnProjections(TracedPoints, "TracedPoints");
		//}
		if (!giDetectSoma)
			gTheVessels.Print(const_cast<char*> (SummaryFileName.c_str()));
		
	}

	//char achTempBuff[128];
	//sprintf(achTempBuff, "%s%s.asc", output_path.c_str(), image_name.c_str());



	//cout << "Neurolucida file: " << achTempBuff << endl; 

	//gLogFile.close();

	FreeResources();

	cout << "Overall Execution Time: " << overall_time.elapsedSec()
		<< " seconds" << endl;
	gLogFile << "Overall Execution Time: " << overall_time.elapsedSec()
		<< " seconds" << endl;
	return 0;
}

