#include "SomaExtraction.h"
#include "ftkUtils.h"
#include <itkNeighborhood.h>
#ifdef _OPENMP
#include "omp.h"
#endif

void MakeDices(const char * montagefileName, const char * seedfileName, int diceWidth,
			   double alfa, double beta, int timethreshold, double curvatureScaling, double rmsThres, int holeSize, int minObjSize);

int main(int argc, char* argv[])
{
	std::string str = std::string("XML");
	if( argc == 5 && std::string(argv[1]) == str && atoi(argv[2]) >= 1)  // generate project file for the image sequences
	{
		SomaExtractor *Somas = new SomaExtractor();	
		SomaExtractor::OutputImageType::Pointer inputImage = Somas->Read8BitImage("1_8bit.tif");
		int width = inputImage->GetLargestPossibleRegion().GetSize()[0];
		int height = inputImage->GetLargestPossibleRegion().GetSize()[1];
		int row = atoi(argv[3]);   // arrange them row by col
		int col = atoi(argv[4]);
		std::ofstream ofs("Project.xml");
		ofs<< "<?xml	version=\"1.0\"	?>"<<std::endl;
		ofs<< "<Source>"<<std::endl;
		for( int i = 1; i <= atoi(argv[2]); i++)
		{
			int rown = (i-1) / col;
			int coln = (i-1) % col;
			int tx = width * 1.1 * coln;
			int ty = height * 1.1 * rown;
			ofs<<"\t"<<"<File	FileName=\""<<i<<"_CV_ANT.swc\""<<"\t"<<"Type=\"Trace\"\ttX=\""<<tx<<"\"\ttY=\""<<ty<<"\"\ttZ=\"0\"/>"<<std::endl;
			ofs<<"\t"<<"<File	FileName=\""<<i<<"_8bit.tif\""<<"\t"<<"Type=\"Image\"\ttX=\""<<tx<<"\"\ttY=\""<<ty<<"\"\ttZ=\"0\"/>"<<std::endl;
			ofs<<"\t"<<"<File	FileName=\""<<i<<"_soma.mhd\""<<"\t"<<"Type=\"Soma\"\ttX=\""<<tx<<"\"\ttY=\""<<ty<<"\"\ttZ=\"0\"/>"<<std::endl;
			//ofs<<"\t"<<"<File	FileName=\""<<i<<"_soma_features.txt\""<<"\t"<<"Type=\"Nuclei_Table\"\ttX=\""<<tx<<"\"\ttY=\""<<ty<<"\"\ttZ=\"0\"/>"<<std::endl;
		}
		ofs<< "</Source>"<<std::endl;
		delete Somas;
		return 0;
	}

	if( argc < 2 || (atoi(argv[1]) < 0 && atoi(argv[1]) > 2))
	{
		std::cout<<"Debris: SomaExtraction <0> <IntensityImage> <DebrisImage> <SomaSeeds.txt>"<<std::endl;
		//std::cout<<"Derbis: SomaExtraction <InputImageFileName> <Centroids.txt> <DiceWidth (typically 100)> <hole filling (typically 10)>\n";
		std::cout<<"SomaExtraction: SomaExtraction <1> <InputImageFileName> <InitialContourLabeledImage> <Options>\n";
		std::cout<<"SomaExtraction without seeds: SomaExtraction <2> <InputImageFileName> <SomaCentroids.txt> <Options>\n";
		std::cout<<"Normalize intensity: SomaExtraction <3> <InputImageFileName> <BackgroundImageFileName>\n";
		return 0;
	}

	SomaExtractor *Somas = new SomaExtractor();	
	if( atoi(argv[1]) == 0)   /// debris accumulation  
	{
		std::cout<< "Reading Montage1"<<std::endl;
		SomaExtractor::OutputImageType::Pointer inputImage = Somas->Read8BitImage(argv[2]);
		std::cout<< "Reading Montage2"<<std::endl;
		SomaExtractor::OutputImageType::Pointer debrisImage = Somas->Read8BitImage(argv[3]);
		std::vector< itk::Index<3> > somaSeeds;
		std::vector< itk::Index<3> > debrisSeeds;

		Somas->ReadSeedpoints( argv[4], somaSeeds, 0);
		std::cout<< somaSeeds.size()<<std::endl;
		std::cout<< "Generating Debris Centroids..."<<std::endl;
		Somas->GetDebrisCentroids( debrisImage, debrisSeeds);
		//Somas->ReadSeedpoints( argv[3], debrisSeeds, 0);
		std::cout<< debrisSeeds.size()<<std::endl;
		std::cout<< "Associating Debris with Nucleus..."<<std::endl;
		Somas->AssociateDebris(inputImage, somaSeeds, debrisSeeds);
	}
	else if( atoi(argv[1]) == 1) /// soma segmentation
	{
		Somas->LoadOptions( argv[4]); // Load params
		
		std::cout << "Set input image" << std::endl;
		//SomaExtractor::ProbImageType::Pointer image = Somas->SetInputImageByPortion(argv[1]); // Load microglia image by portion
		SomaExtractor::ProbImageType::Pointer image = Somas->SetInputImage(argv[2]); // Load microglia image 16bit

		std::string InputFilename = std::string(argv[2]);
		std::string somaImageName = InputFilename;
		somaImageName.erase(somaImageName.length()-4,somaImageName.length());
		somaImageName.append("_soma.mhd");
		std::string somaCentroidName = InputFilename;
		somaCentroidName.erase(somaCentroidName.length()-4,somaCentroidName.length());
		somaCentroidName.append("_centroids.txt");
		std::string somaFeatureName = InputFilename;
		somaFeatureName.erase(somaFeatureName.length()-4,somaFeatureName.length());
		somaFeatureName.append("_soma_features.txt");

		std::cout << "Set initial contour" << std::endl;
		//SomaExtractor::SegmentedImageType::Pointer initialContourImage = Somas->SetInitalContourImageByPortion(argv[2]); // Load labeled nucleus image by portion 16 bit
		SomaExtractor::SegmentedImageType::Pointer initialContourImage = Somas->SetInitalContourImage(argv[3]); // Load labeled nucleus image
		
		std::vector< itk::Index<3> > seedVector;
		//image = Somas->EnhanceContrast(image, seedVector[0][2], atof(argv[4]), atof(argv[5]));

		clock_t SomaExtraction_start_time = clock();
		
		std::cout<< "Segmenting..."<<std::endl;

		/// SegmentSoma2: GVF Active Contour
		SomaExtractor::SegmentedImageType::Pointer segImage = Somas->SegmentSoma(image, initialContourImage, seedVector);
		std::cout << "Total time for SomaExtraction is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

		/// Compute soma features and write new seeds back
		if( segImage)
		{
			std::cout<< "Writing "<< somaImageName<<std::endl;
			Somas->writeImage(somaImageName.c_str(), segImage);
			std::cout<< "Writing "<< somaCentroidName<<std::endl;
			Somas->writeCentroids( somaCentroidName.c_str() ,seedVector);
			vtkSmartPointer<vtkTable> table = Somas->ComputeSomaFeatures(segImage);
			std::cout<< "Writing "<< somaFeatureName<<std::endl;
			ftk::SaveTable(somaFeatureName.c_str(), table);
		}
	}
	else if( atoi(argv[1]) == 2)  /// soma segmentation without initial seeds
	{
		Somas->LoadOptions( argv[4]); // Load params
		std::cout << "Set input image" << std::endl;
		SomaExtractor::OutputImageType::Pointer image = Somas->Read16BitImage(argv[2]); // Load microglia image 16bit	
		
		std::string InputFilename = std::string(argv[2]);
		std::string bit8FileName = InputFilename;
		bit8FileName.erase(bit8FileName.length()-4,bit8FileName.length());
		bit8FileName.append("_8bit.tif");
		Somas->writeImage(bit8FileName.c_str(), image);

		std::string somaImageName = InputFilename;
		somaImageName.erase(somaImageName.length()-4,somaImageName.length());
		somaImageName.append("_soma.mhd");
		std::string somaCentroidName = InputFilename;
		somaCentroidName.erase(somaCentroidName.length()-4,somaCentroidName.length());
		somaCentroidName.append("_centroids.txt");
		std::string somaFeatureName = InputFilename;
		somaFeatureName.erase(somaFeatureName.length()-4,somaFeatureName.length());
		somaFeatureName.append("_soma_features.txt");

		std::vector< itk::Index<3> > seedVector;
		SomaExtractor::ProbImageType::Pointer binImagePtr = Somas->GenerateSeedPoints(image, seedVector);
		Somas->ReadSeedpoints(argv[3], seedVector, false);

		clock_t SomaExtraction_start_time = clock();
		
		std::cout<< "Segmenting..."<<std::endl;

		/// SegmentSoma1: Active Contour without GVF, eliminate small objects
		SomaExtractor::SegmentedImageType::Pointer segImage = Somas->SegmentSoma(image, seedVector, binImagePtr);
		std::cout << "Total time for SomaExtraction is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

		/// Compute soma features and write new seeds back
		if( segImage)
		{
			std::cout<< "Writing "<< somaImageName<<std::endl;
			Somas->writeImage(somaImageName.c_str(), segImage);
			std::cout<< "Writing "<< somaCentroidName<<std::endl;
			Somas->writeCentroids( somaCentroidName.c_str() ,seedVector);
			vtkSmartPointer<vtkTable> table = Somas->ComputeSomaFeatures(segImage);
			std::cout<< "Writing "<< somaFeatureName<<std::endl;
			ftk::SaveTable(somaFeatureName.c_str(), table);
		}
	}
	else if( atoi(argv[1]) == 3 && argc >= 4 && argc <= 5)  /// normalize the intensity
	{
		std::string InputFilename = std::string(argv[2]);
		std::string imageName = InputFilename;
		imageName.erase(imageName.length()-4,imageName.length());
		imageName.append("_normalize.tif");

		SomaExtractor::ProbImageType::Pointer image = Somas->SetInputImage(argv[2]); // Load microglia image 16bit
		SomaExtractor::ProbImageType2D::Pointer imageBackground = Somas->SetInputImage2D(argv[3]); // Load microglia image 16bit
		// Devide image by imageBackground, rescale the result to the maximum scale of the origional image scale
		SomaExtractor::UShortImageType::Pointer resultImage;
		if( argc == 5 && atoi(argv[4]) > 1)
		{
			resultImage = Somas->DevideAndScale(image, imageBackground, atoi(argv[4])); 
		}
		else
		{
			resultImage = Somas->DevideAndScale(image, imageBackground); 
		}
		Somas->writeImage(imageName.c_str(), resultImage);
	}
	delete Somas;
	return 0;
}

void MakeDices(const char * montagefileName, const char * seedfileName, int diceWidth, 
			   double alfa, double beta, int timethreshold, double curvatureScaling, double rmsThres, int holeSize, int minObjSize)
{
	typedef SomaExtractor::ProbImageType ImageType;
	typedef SomaExtractor::OutputImageType OutputImageType;
	typedef itk::RegionOfInterestImageFilter< ImageType, ImageType> RegionOfInterestFilter;
	SomaExtractor *Somas = new SomaExtractor();

	std::cout<< "ReadMontage..."<<std::endl;
	ImageType::Pointer image = Somas->SetInputImage(montagefileName);
	int SZ = image->GetLargestPossibleRegion().GetSize()[2];

	std::vector< itk::Index<3> > seedVector;
	Somas->ReadSeedpoints( seedfileName, seedVector, 0);
	std::cout<< "Seed points size:"<< seedVector.size()<<std::endl;	

	#pragma omp parallel for
	for(int i = 0; i < seedVector.size(); i++)
	{
		SomaExtractor *somaExtractor = new SomaExtractor();
		ImageType::IndexType start;
		start[0] = seedVector[i][0] - diceWidth / 2;
		start[1] = seedVector[i][1] - diceWidth / 2;
		start[2] = 0;
		ImageType::SizeType size;
		size[0] = diceWidth;
		size[1] = diceWidth;
		size[2] = SZ;

		ImageType::RegionType desiredRegion;
		desiredRegion.SetSize(size);
		desiredRegion.SetIndex(start);

		RegionOfInterestFilter::Pointer regionFilter = RegionOfInterestFilter::New();
		regionFilter->SetInput(image);
		regionFilter->SetRegionOfInterest(desiredRegion);
		regionFilter->Update();
		ImageType::Pointer diceImage = regionFilter->GetOutput();
		std::ostringstream oss1;
		oss1<< "Dice_"<< i<<".tif";
		std::string str1 = oss1.str();
		somaExtractor->writeImage(str1.c_str(), diceImage);
		
		std::vector< itk::Index<3> > seedsForDice;
		itk::Index<3> seedIndex = seedVector[i];
		seedIndex[0] = diceWidth / 2;
		seedIndex[1] = diceWidth / 2;
		seedsForDice.push_back(seedIndex);

		float Origion[3] = {0,0,0};
		diceImage->SetOrigin(Origion);

		//SomaExtractor::SegmentedImageType::Pointer segImage= somaExtractor->SegmentSoma(diceImage, seedsForDice, alfa, 
		//					beta, timethreshold, curvatureScaling, rmsThres, holeSize, minObjSize);
		//double threshold = 0;
//		//ImageType::Pointer segImage = Somas->EnhanceContrast(diceImage, seedIndex[2], alfa, beta, threshold);
//
//#pragma omp critical
//		std::cout<< "segment "<< i<<"\t"<<seedsForDice.size()<<std::endl;
//		std::ostringstream oss;
//		oss<< "Dice_Segment_"<< i<<".tif";
//		std::string str = oss.str();
//		somaExtractor->writeImage(str.c_str(), segImage);
		delete somaExtractor;
	}
	delete Somas;
}