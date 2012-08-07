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
	if( argc != 4 && argc != 7)
	{
		std::cout<<"Debris: SomaExtraction <IntensityImage> <DebrisImage> <SomaSeeds.txt>"<<std::endl;
		//std::cout<<"Derbis: SomaExtraction <InputImageFileName> <Centroids.txt> <DiceWidth (typically 100)> <hole filling (typically 10)>\n";
		std::cout<<"SomaExtraction: SomaExtraction <InputImageFileName> <InitialContourLabeledImage> <SomaCentroids.txt> <Nucleus Table:1; Centroids table:0> <Options> <SomaImage.tif> \n";
		return 0;
	}

	/// debris accumulation
	if( argc == 4)    
	{
		SomaExtractor *Somas = new SomaExtractor();
		std::cout<< "Reading Montage1"<<std::endl;
		SomaExtractor::OutputImageType::Pointer inputImage = Somas->Read8BitImage(argv[1]);
		std::cout<< "Reading Montage2"<<std::endl;
		SomaExtractor::OutputImageType::Pointer debrisImage = Somas->Read8BitImage(argv[2]);
		std::vector< itk::Index<3> > somaSeeds;
		std::vector< itk::Index<3> > debrisSeeds;

		Somas->ReadSeedpoints( argv[3], somaSeeds, 0);
		std::cout<< somaSeeds.size()<<std::endl;
		std::cout<< "Generating Debris Centroids..."<<std::endl;
		Somas->GetDebrisCentroids( debrisImage, debrisSeeds);
		//Somas->ReadSeedpoints( argv[3], debrisSeeds, 0);
		std::cout<< debrisSeeds.size()<<std::endl;
		std::cout<< "Associating Debris with Nucleus..."<<std::endl;
		Somas->AssociateDebris(inputImage, somaSeeds, debrisSeeds);
		delete Somas;
	}

	/// soma segmentation
	if( argc == 7)
	{
		SomaExtractor *Somas = new SomaExtractor();

		std::cout << "Entering SetInputImage" << std::endl;
		SomaExtractor::ProbImageType::Pointer image = Somas->SetInputImage(argv[1]); // Load microglia image
		SomaExtractor::SegmentedImageType::Pointer initialContourImage = Somas->SetInitalContourImage(argv[2]); // Load labeled nucleus image
		std::vector< itk::Index<3> > seedVector;
		Somas->LoadOptions( argv[5]); // Load params

		bool bNucleusTable = false;
		if(atoi(argv[4]) != 0)
		{
			bNucleusTable = true;
		}
		
		Somas->ReadSeedpoints( argv[3], seedVector, bNucleusTable);

		//image = Somas->EnhanceContrast(image, seedVector[0][2], atof(argv[4]), atof(argv[5]));

		clock_t SomaExtraction_start_time = clock();
		
		std::cout<< "Segmenting..."<<std::endl;

		/// SegmentSoma2: GVF Active Contour
		SomaExtractor::SegmentedImageType::Pointer segImage = Somas->SegmentSoma2(image, initialContourImage, seedVector);
		std::cout << "Total time for SomaExtraction is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

		/// Compute soma features and write new seeds back
		if( segImage)
		{
			std::cout<< "Writing Soma Image."<<std::endl;
			Somas->writeImage(argv[6], segImage);
			Somas->writeCentroids( "NewSomaCentroids.txt" ,seedVector);

			vtkSmartPointer<vtkTable> table = Somas->ComputeSomaFeatures(segImage);
			ftk::SaveTable("SomaFeatures.txt", table);

			//std::cout<< "Writing Soma Seeds Image."<<std::endl;
			//Somas->WriteSomaSeedsIntoImage();

			std::cout<< "Writing Soma Features."<<std::endl;
			if(bNucleusTable)
			{
				ftk::SaveTable("SomaFeaturesWithAllSeeds.txt", table);
			}
			else
			{
				ftk::SaveTable("SomaFeatures.txt", table);
			}
		}
		delete Somas;
	}
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