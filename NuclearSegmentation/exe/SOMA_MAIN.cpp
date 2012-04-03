#include "SomaExtraction.h"
#include "ftkUtils.h"
#include <itkNeighborhood.h>
#ifdef _OPENMP
#include "omp.h"
#endif

void MakeDices(const char * montagefileName, const char * seedfileName, const char * outputDir, int diceWidth,
			   double alfa, double beta, int timethreshold, double curvatureScaling, double rmsThres, int minObjSize);

int main(int argc, char* argv[])
{	
	if( argc != 5 && argc != 11)
	{
		std::cout<<"MakeDices: SomaExtraction <InputImageFileName> <Centroids.txt> <OutputDir> <DiceWidth (typically 100)>\n";
		std::cout<<"SomaExtraction: SomaExtraction <InputImageFileName> <SomaCentroids.txt> <SomaImage.tif> <alfa 1> <beta 30> <Time Threshold (typically near 10)> <Curvature Scaling (typically 0.5)>\
					<RMS Error (typically 0.02)> <Min Object Size (typically 1500)> <Nucleus Table:1; Centroids table:0>\n";
		return 0;
	}

	if( argc == 5)
	{
		MakeDices(argv[1], argv[2], argv[3], atoi(argv[4]), 2.0, 35.0, 5, 0.6, 0.02, 1000);
	}

	if( argc == 11)
	{
		SomaExtractor *Somas = new SomaExtractor();

		std::cout << "Entering SetInputImage" << std::endl;
		Somas->SetInputImage(argv[1]);
		std::vector< itk::Index<3> > seedVector;

		bool bNucleusTable = false;
		if(atoi(argv[10]) != 0)
		{
			bNucleusTable = true;
		}

		Somas->ReadSeedpoints( argv[2], seedVector, bNucleusTable);
		//SomaExtractor::ProbImageType::Pointer image = Somas->GetEdgePotentialMap();
		SomaExtractor::ProbImageType::Pointer image = Somas->GetFloatInputImage();

		std::cout<< seedVector.size()<<std::endl;

		image = Somas->EnhanceContrast(image, atof(argv[4]), atof(argv[5]));

		clock_t SomaExtraction_start_time = clock();
		SomaExtractor::SegmentedImageType::Pointer segImage = Somas->SegmentSoma(image, seedVector, atof(argv[4]), atof(argv[5]), atoi(argv[6]), atof(argv[7]), atof(argv[8]), atoi(argv[9]));

		std::cout<< "Writing Soma Image."<<std::endl;
		Somas->writeImage(argv[3], segImage);
		Somas->writeCentroids( "NewSomaCentroids.txt" ,seedVector);

		vtkSmartPointer<vtkTable> table = Somas->ComputeSomaFeatures(segImage);
		ftk::SaveTable("SomaFeatures.txt", table);

		////std::cout<< "Writing Soma Seeds Image."<<std::endl;
		////Somas->WriteSomaSeedsIntoImage();

		//vtkSmartPointer<vtkTable> table = Somas->GetSomaFeatureTable();

		//std::cout<< "Writing Soma Features."<<std::endl;
		//if(bNucleusTable)
		//{
		//	ftk::SaveTable("SomaFeaturesWithAllSeeds.txt", table);
		//}
		//else
		//{
		//	ftk::SaveTable("SomaFeatures.txt", table);
		//}
		std::cout << "Total time for SomaExtraction is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;
		delete Somas;
	}
	return 0;
}

void MakeDices(const char * montagefileName, const char * seedfileName, const char * outputDir, int diceWidth, 
			   double alfa, double beta, int timethreshold, double curvatureScaling, double rmsThres, int minObjSize)
{
	typedef SomaExtractor::ProbImageType ImageType;
	typedef SomaExtractor::OutputImageType OutputImageType;
	typedef itk::RegionOfInterestImageFilter< ImageType, ImageType> RegionOfInterestFilter;
	SomaExtractor *Somas = new SomaExtractor();

	std::cout<< "ReadMontage..."<<std::endl;
	Somas->SetInputImage(montagefileName);
	ImageType::Pointer image = Somas->GetFloatInputImage();
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

		std::vector< itk::Index<3> > seedsForDice;
		itk::Index<3> seedIndex = seedVector[i];
		seedIndex[0] = diceWidth / 2;
		seedIndex[1] = diceWidth / 2;
		seedsForDice.push_back(seedIndex);
		std::cout<< "Segment "<< i<<std::endl;

		float Origion[3] = {0,0,0};
		diceImage->SetOrigin(Origion);

		SomaExtractor::SegmentedImageType::Pointer segImage= somaExtractor->SegmentSoma(diceImage, seedsForDice, alfa, 
							beta, timethreshold, curvatureScaling, rmsThres, minObjSize);

		std::ostringstream oss;
		oss<< "Dice_Seg_"<< i<<".tif";
		std::string str = oss.str();
		somaExtractor->writeImage(str.c_str(), segImage);
		delete somaExtractor;
	}
	delete Somas;
}