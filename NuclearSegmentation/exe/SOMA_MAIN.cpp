#include "SomaExtraction.h"
#include "ftkUtils.h"

int main(int argc, char* argv[])
{	
	clock_t SomaExtraction_start_time = clock();

	if(argc < 7)
	{
		std::cout<<"Usage: SomaExtraction <InputImageFileName> <SomaCentroids.txt> <SomaImage.tif>\
			<Time Threshold (typically near 30)> <Curvature Scaling (typically 0.5)> <RMS Error (typically 0.02)>\n";
		return 0;
	}

	SomaExtractor *Somas = new SomaExtractor();
	
	std::cout << "Entering SetInputImage" << std::endl;
	Somas->SetInputImage(argv[1]);
	std::vector< itk::Index<3> > seedVector;
	Somas->ReadSeedpoints( argv[2], seedVector);
	
	SomaExtractor::SegmentedImageType::Pointer segImage = Somas->SegmentSoma(Somas->GetFloatInputImage(), seedVector, atoi(argv[4]), atof(argv[5]), atof(argv[6]));

	Somas->writeSomaImage(argv[3]);

	vtkSmartPointer<vtkTable> table = Somas->ComputeSomaFeatures(segImage);
	ftk::SaveTable("SomaFeatures.txt", table);

	std::cout << "Total time for SomaExtraction is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

	return 0;
}