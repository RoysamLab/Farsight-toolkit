#include "SomaExtraction.h"
#include "ftkUtils.h"

int main(int argc, char* argv[])
{	
	clock_t SomaExtraction_start_time = clock();

	if(argc < 10)
	{
		std::cout<<"Usage: SomaExtraction <InputImageFileName> <SomaCentroids.txt> <SomaImage.tif> <Time Threshold (typically near 10)> <Curvature Scaling (typically 0.5)>\
				    <RMS Error (typically 0.02)> <Min Object Size (typically 1500)> <Volume Threshold (typically 5000)> <Nucleus Table:1; Centroids table:0>\n";
		return 0;
	}

	SomaExtractor *Somas = new SomaExtractor();
	
	std::cout << "Entering SetInputImage" << std::endl;
	Somas->SetInputImage(argv[1]);
	std::vector< itk::Index<3> > seedVector;

	bool bNucleusTable = false;
	if(atoi(argv[9]) != 0)
	{
		bNucleusTable = true;
	}

	Somas->ReadSeedpoints( argv[2], seedVector, bNucleusTable);
	//Somas->SmoothByDiffusionImageFilter();
	
	SomaExtractor::SegmentedImageType::Pointer segImage = Somas->SegmentSoma(Somas->GetFloatInputImage(), seedVector, atoi(argv[4]), atof(argv[5]), atof(argv[6]), atoi(argv[7]), atoi(argv[8]), bNucleusTable);

	std::cout<< "Writing Soma Image."<<std::endl;
	Somas->writeSomaImage(argv[3]);

	//std::cout<< "Writing Soma Seeds Image."<<std::endl;
	//Somas->WriteSomaSeedsIntoImage();

	vtkSmartPointer<vtkTable> table = Somas->GetSomaFeatureTable();

	std::cout<< "Writing Soma Features."<<std::endl;
	if(bNucleusTable)
	{
		ftk::SaveTable("SomaFeaturesWithAllSeeds.txt", table);
	}
	else
	{
		ftk::SaveTable("SomaFeatures.txt", table);
	}

	std::cout << "Total time for SomaExtraction is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

	return 0;
}