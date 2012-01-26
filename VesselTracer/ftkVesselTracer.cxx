
#include "ftkVesselTracer.h"

ftkVesselTracer::ftkVesselTracer(){
}

ftkVesselTracer::ftkVesselTracer(std::string input_data_path, bool preprocess = false, bool start_with_mst = false){

	if(preprocess){
		this->allParams.preProcessingParams.initByDefaultValues();
		ImageType3D::Pointer tiff_data_ptr; //only required for preprocessing
		this->PreprocessData(input_data_path, tiff_data_ptr); // ask user	
	}

	// Load preprocessed data
	this->LoadPreprocessedData(input_data_path);

	if(start_with_mst == false){

		// Computing spheical bins. This could be a part of preprocessing - need to write the output structure in a file.
		this->allParams.oriBin.initByDefaultValues();
		this->SphericalBinPreprocess();

		// 3D renering of original and preprocessed data
		Common::RescaleDataForRendering(this->originalData, this->originalDataForRendering);
		Common::RescaleDataForRendering(this->inputData, this->inputDataForRendering);		
		//Common::RenderImage3D(this->originalDataForRendering);
		//Common::RenderImage3D(this->inputDataForRendering);
		

		// MIP computation
		this->ComputeIntensityProjectionImages();
		//this->RenderMaximumProjectionImage();	

		// Primary and secondary node detection
		this->allParams.nodeDetectionParams.initByDefaultValues();
		this->ComputeAllPrimaryNodes();
		this->ComputeAllSecondaryNodes();

		// MST and post processing
		this->allParams.graphAndMSTParams.initByDefaultValues();
		this->CreateMinimumSpanningForest();
	}
	else{

		Common::RescaleDataForRendering(this->originalData, this->originalDataForRendering);
		Common::RescaleDataForRendering(this->inputData, this->inputDataForRendering);
		
		// Reading secondary nodes from file for now
		std::string filename = "AllNodes_grid10.txt";
		this->ReadNodesFromTextFile(filename);
		//this->ReadNodesFromTextFile(std::string("AllNodes_default.txt"));
		
		this->allParams.graphAndMSTParams.initByDefaultValues();
		this->CreateMinimumSpanningForest();
	}
}
ftkVesselTracer::~ftkVesselTracer(){
}

void PreprocessingParameters::initByDefaultValues(void){
	
	this->medianFilterRadius = 1;
	this->anisDiffusionNIter = 10;
	this->anisDiffusionConductance = 30;
	this->smoothingSigma = 2.0f;
}

void SphericalBinInfo::initByDefaultValues(void){
	
	this->indexLength = 50;
	this->angleCount = 182;
	this->angleInrement = 10;
	this->nLastIndicesOfInterest = 8;
	this->histSmoothingFactor = 0.05;
	this->minSphHistCount = 0.0001;
}

void NodeDetectionParameters::initByDefaultValues(void){
	
	this->gridSpacing = 30; //10; //20; //30; //10; //15; //20; //15;
	this->iterNPrimaryNode = 100;

	this->increaseLikelihoodThreshold = 0.001;
	this->discardNodeLikelihoodThreshold = 0.01; // 0.001 accroding to thesis
	this->iterNForOnlyRegionBasedTerm = 15;
	this->iterNMinimum = 10;
	this->regionBasedTermWeight = 0.3; // 30%
	this->edgeBasedTermWeight = 0.7; // 70%
	this->minimumAccumulatedParamChange = 0.1;
	this->iterNMonitorParamChange = 5;
	this->currentMonitoredIter = 0;

	this->dtX = 5000; //50000; 
	this->dtY = 5000; //50000; 
	this->dtZ = 5000; //50000;
	this->dtScale = 1000; //10000;
	this->dirX = 0; this->dirY = 0; this->dirZ = 0;
	this->dirScale = 0;
	this->chX = std::vector<double>(this->iterNMonitorParamChange, 0); 
	this->chY = std::vector<double>(this->iterNMonitorParamChange, 0); 
	this->chZ = std::vector<double>(this->iterNMonitorParamChange, 0);
	this->chScale = std::vector<double>(this->iterNMonitorParamChange, 0);
	
	this->maxVesselWidth = 20;
	this->minVesselWidth = 2;
	this->likelihoodThresholdPrimary = 0.01; //0.005;
	this->distanceThresholdPrimary = 1.2;

	this->traceQualityThreshold = 4.0;
	this->maxQueueIter = 20000;
	this->distanceThresholdSecondary = 1.0;
	this->maxBranchAngle = 60; // in degrees
	this->branchingThreshold = 0.05;
	this->secondaryParamRateReduction = 1000;
	this->dtXSecondary = this->dtX / this->secondaryParamRateReduction;
	this->dtYSecondary = this->dtY / this->secondaryParamRateReduction;
	this->dtZSecondary = this->dtZ / this->secondaryParamRateReduction;
	this->dtScaleSecondary = this->dtScale / this->secondaryParamRateReduction;
	this->iterNMinimumSecondary = 15;
	this->regionBasedTermWeightSecondary = 0.5;
	this->edgeBasedTermWeightSecondary = 0.5;
	this->secondarySearchConstraint = std::sqrt(2.0);
	this->primaryReversePositionRate = 0.5;
	this->primaryReverseScaleRate = 0.5;
	this->secondaryReversePositionRate = 0.8;
	this->secondaryReverseScaleRate = 0.9;
	this->maxTraceCost = 10.0;
	this->infTraceQuality = 100;
	this->maxQueueSize = 10000;
	this->traceLengthCost = 0.5; //1.0; //0.5; //1.0; //0.5; //0.0;
	this->primaryNodeSearchRadFactor = 0.75;
}

void GraphAndMSTPartameters::initByDefaultValues(void){

	this->affinityRadThresh = 2.0;
	this->NBinsAffinity = 32;
	this->maxEdgeWeight = 99.0;
	this->minBranchAngle = vnl_math::pi/8;
	this->maxNBranches = 25; //10; //5;
	this->maxTreeNodes = 100;
}
void AllParameters::initByDefaultValues(void){
	
	this->preProcessingParams.initByDefaultValues();
	this->oriBin.initByDefaultValues();
	this->nodeDetectionParams.initByDefaultValues();
	this->graphAndMSTParams.initByDefaultValues();
}

int ftkVesselTracer::PreprocessData(std::string file_path, ImageType3D::Pointer& data_ptr){

	itk::TimeProbe timer; // Timer for measuring the preprocessing time
	timer.Start();
	
	std::string write_file_path = file_path.substr(0, file_path.find_last_of('.'));
	std::string append_str = "_original.mhd";

	try{
		Common::ReadImage3D(file_path, data_ptr);

		Common::WriteImage3D(write_file_path + std::string(append_str), data_ptr);
	}
	catch(itk::ExceptionObject& e){
		std::cout << e << std::endl;
		return EXIT_FAILURE;
	}

	//PARAMS: Make available outside function
	int median_radius = this->allParams.preProcessingParams.medianFilterRadius;
	int anis_diffusion_N_iter = this->allParams.preProcessingParams.anisDiffusionNIter;
	int anis_diffusion_conductance = this->allParams.preProcessingParams.anisDiffusionConductance;
	float smoothing_sigma = this->allParams.preProcessingParams.smoothingSigma;

	Common::MedianFilter(median_radius, data_ptr);
	Common::CurvatureAnisotropicDiffusion(anis_diffusion_N_iter, anis_diffusion_conductance, data_ptr);
	append_str = "_preprocessed.tif";

	try{
		Common::WriteTIFFImage3D(write_file_path + std::string(append_str), data_ptr);
	}
	catch(itk::ExceptionObject& e){
		std::cout << e << std::endl;
		return EXIT_FAILURE;
	}

	Common::GVFDiffusion(smoothing_sigma, write_file_path, data_ptr);	

	// Compute oriBin and save to file?

	timer.Stop();
	std::cout << "The processing took " << timer.GetMeanTime() << " seconds. " << std::endl;

	return EXIT_SUCCESS;
}

void ftkVesselTracer::LoadPreprocessedData(std::string data_path){
	
	data_path = data_path.substr(0, data_path.find_last_of('.')); // get rid of the extension
	std::string append_ext = ".mhd";
	std::string append_original = "_original.mhd";
	std::string append_gx = "_gx.mhd";
	std::string append_gy = "_gy.mhd";
	std::string append_gz = "_gz.mhd";

	try{
		Common::ReadImage3D(data_path + std::string(append_ext), this->inputData); //input data is preprocessed data

		Common::ReadImage3D(data_path + std::string(append_original), this->originalData);	
		Common::ReadImage3D(data_path + std::string(append_gx), this->gx);
		Common::ReadImage3D(data_path + std::string(append_gy), this->gy);
		Common::ReadImage3D(data_path + std::string(append_gz), this->gz);
	}
	catch(itk::ExceptionObject& e){
		std::cout << e << std::endl;		
	}

	//testing the data..
	/*MinMaxCalculatorType::Pointer min_max_filter = MinMaxCalculatorType::New();
	//min_max_filter->SetImage(this->inputData);
	//min_max_filter->SetImage(this->originalData);
	min_max_filter->SetImage(this->gx);
	min_max_filter->Compute();

	float min = min_max_filter->GetMinimum();
	float max = min_max_filter->GetMaximum();*/

	/*InvertImageFilterType::Pointer image_inverter = InvertImageFilterType::New();
	image_inverter->SetInput(this->originalData);
	image_inverter->Update();
	
	this->originalData = image_inverter->GetOutput();*/
	
	
	// Same pixel apparently has different values in ITK and Matlab!!
	/*ImageType3D::IndexType test_index, test_index1;
	test_index[0] = 123; test_index[1] = 6; test_index[2] = 23;
	test_index1[0] = 6; test_index1[1] = 123; test_index1[2] = 23;
	PixelType test_pixel = this->inputData->GetPixel(test_index);
	PixelType test_pixel1 = this->inputData->GetPixel(test_index1);*/
		
}

void ftkVesselTracer::ComputeIntensityProjectionImages(void){
	
	MaxProjectionFilterType::Pointer max_intensity_projector = MaxProjectionFilterType::New();
	max_intensity_projector->SetInput(this->inputDataForRendering);

	max_intensity_projector->Update();

	this->maximumProjectionImage = max_intensity_projector->GetOutput();

	MinProjectionFilterType::Pointer min_intensity_projector = MinProjectionFilterType::New();
	min_intensity_projector->SetInput(this->inputDataForRendering);

	min_intensity_projector->Update();

	this->minimumProjectionImage = min_intensity_projector->GetOutput();


	//CANNOT RETURN VTK IMAGE AFTER USING ITKTOVTK FILTER
	/*RenderImageType3D::Pointer max_intensity_image = max_intensity_projector->GetOutput();

	ITKToVTKConnectorType::Pointer itk_to_vtk_connector = ITKToVTKConnectorType::New();
	
	//itk_to_vtk_connector->SetInput(max_intensity_image);
	itk_to_vtk_connector->SetInput(this->inputDataForRendering);
	
	itk_to_vtk_connector->Update();

	//vtk_image = vtkSmartPointer<vtkImageData>::New();
	//this->MIP_image = itk_to_vtk_connector->GetOutput();

	//vtkSmartPointer<vtkImageActor> vtk_image_actor = vtkSmartPointer<vtkImageActor>::New();

	vtkSmartPointer<vtkImageData> vtk_image = itk_to_vtk_connector->GetOutput();

	MIP_image_actor = vtkSmartPointer<vtkImageActor>::New();
	
	MIP_image_actor->SetInput(vtk_image);
	MIP_image_actor->Update();

	//vtk_image->UpdateInformation();

	//testing projection image
	//vtk_image->PrintSelf(std::cout, vtkIndent(0));

	vtkSmartPointer<vtkImageViewer2> image_viewer = vtkSmartPointer<vtkImageViewer2>::New();
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	image_viewer->SetInput(MIP_image_actor->GetInput()); // CODE BREAKS HERE
	image_viewer->SetupInteractor(render_window_interactor);
	image_viewer->Render();
	image_viewer->GetRenderer()->ResetCamera();
	render_window_interactor->Start();

	return vtk_image;*/
}

void ftkVesselTracer::RenderMaximumProjectionImage(void){

	ITKToVTKConnectorType::Pointer itk_to_vtk_connector = ITKToVTKConnectorType::New();
	
	//itk_to_vtk_connector->SetInput(max_intensity_image);
	itk_to_vtk_connector->SetInput(this->maximumProjectionImage);
	
	itk_to_vtk_connector->Update();

	// The new rendering pipeline
	vtkSmartPointer<vtkImageViewer2> image_viewer = vtkSmartPointer<vtkImageViewer2>::New();
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	
	image_viewer->SetInput(itk_to_vtk_connector->GetOutput()); 
	image_viewer->SetupInteractor(render_window_interactor);
	image_viewer->Render();
	image_viewer->GetRenderer()->ResetCamera();
	render_window_interactor->Start();
	
	// The traditional rendering pipeline
	/*vtkSmartPointer<vtkImageActor> image_actor = vtkSmartPointer<vtkImageActor>::New();
	
	image_actor->SetInput(vtk_image);
	 
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(image_actor);
	renderer->SetBackground(1.0, 1.0, 1.0);
	renderer->ResetCamera();
	
	vtkSmartPointer<vtkRenderWindow> render_window = vtkSmartPointer<vtkRenderWindow>::New();
	render_window->AddRenderer(renderer);
	render_window->SetSize(500, 500);

	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	render_window_interactor->SetRenderWindow(render_window);
	
	render_window->Render();
	
	render_window_interactor->Initialize();
	render_window->Render();

	render_window_interactor->Start();	
	*/
}

bool CompareArrayElementsPreserveIndex2(arrayElement a1, arrayElement a2){
	return(a1.getElement() < a2.getElement());
}

void ftkVesselTracer::SphericalBinPreprocess(void){

	std::cout << "Computing SphericalBinInfo..." << std::endl;

	//ArrayType3D temp(boost::extents[2*indexLength+1][2*indexLength+1][2*indexLength+1]);
	//this->allParams.oriBin.BinIndex = temp; //boost::extents[2*indexLength+1][2*indexLength+1][2*indexLength+1];

	int indexLength = this->allParams.oriBin.indexLength;
	int angleCount = this->allParams.oriBin.angleCount;

	int dim_size = 2 * indexLength + 1;
	
	//VectorType3D temp(dim_size, VectorType2D(dim_size, VectorType1D(dim_size, 0)));
	
	this->allParams.oriBin.binIndexVec = VectorType3D(dim_size, VectorType2D(dim_size, VectorType1D(dim_size, 0)));
	
	std::vector<double> theta, phi, x(angleCount, 0), y(angleCount, 0), z(angleCount, -1);
	
	for(int i = -4; i <= 4; i++)
		theta.push_back(vnl_math::pi * (double)i / 10.0);
	for(int i = -9; i <= 10; i++)	
		phi.push_back(vnl_math::pi * (double)i / 10.0);
	
	//Calculating bin centers
	int count = 1;
	
	//#pragma omp parallel for firstprivate(count)
	for(int i = 0; i < theta.size(); i++){
		for(int j = 0; j < phi.size(); j++){
			x[count] = cos(double(theta[i])) * cos(double(phi[j]));
			y[count] = cos(double(theta[i])) * sin(double(phi[j]));
			z[count] = sin(double(theta[i]));

			count++;
		}
	}
	z[count] = 1;

	// A small test array
	/*double test_point[3] = {0.0, 0.0, 0.0}, test_norm = 1.0;
	int test_max_d_position = 0;
	double test_max_d = 0.0;
	std::vector<double> test_d(x.size(), 1);
	VectorType3D test_vec = VectorType3D(11, VectorType2D(11, VectorType1D(11, 0)));
	for(int i1 = -5; i1 <=5; i1++){
		for(int i2 = -5; i2 <= 5; i2++){
			for(int i3 = -5; i3 <=5; i3++){

				test_point[0] = (double)i1; test_point[1] = (double)i2; test_point[2] = (double)i3;
				test_norm = sqrt((test_point[0]*test_point[0]) + (test_point[1]*test_point[1]) + (test_point[2]*test_point[2]));
				test_point[0] = test_point[0] / test_norm; test_point[1] = test_point[1] / test_norm; test_point[2] = test_point[2] / test_norm;
				
				// Trying to find which point in the cartesian space is closest to the current bin point
				for(int i = 0; i < x.size(); i++)
					test_d[i] = (test_point[0] * x[i]) + (test_point[1] * y[i]) + (test_point[2] * z[i]);

				test_max_d = *std::max_element(test_d.begin(), test_d.end());
				test_max_d_position = std::find(test_d.begin(), test_d.end(), test_max_d) - test_d.begin(); 
				
				test_vec[i3 + 5][i2 + 5][i1 + 5] = test_max_d_position;
			}
		}
	}
	int t1 = test_vec[3][5][6];
	int t2 = test_vec[5][2][9];
	*/

	#pragma omp parallel for 
	for(int i1 = -indexLength; i1 <= indexLength; i1++){
		for(int i2 = -indexLength; i2 <= indexLength; i2++){
			for(int i3 = -indexLength; i3 <= indexLength; i3++){
				
				double a_point[3] = {0.0, 0.0, 0.0}, norm = 1.0;
				a_point[0] = (double)i1; a_point[1] = (double)i2; a_point[2] = (double)i3;
				norm = sqrt((a_point[0]*a_point[0]) + (a_point[1]*a_point[1]) + (a_point[2]*a_point[2]));
				a_point[0] = a_point[0] / norm; a_point[1] = a_point[1] / norm; a_point[2] = a_point[2] / norm;
				
				// Trying to find which point in the cartesian space is closest to the current bin point
				std::vector<double> d(x.size(), 1);
				for(int i = 0; i < x.size(); i++)
					d[i] = (a_point[0] * x[i]) + (a_point[1] * y[i]) + (a_point[2] * z[i]);

				double max_d = *std::max_element(d.begin(), d.end());
				int max_d_position = std::find(d.begin(), d.end(), max_d) - d.begin(); 
				
				// CHECK THE VALUES HERE 
				
				//this->allParams.oriBin.binIndexVec[i1+indexLength+1][i2+indexLength+1][i3+indexLength+1] = max_d_position;
				
				this->allParams.oriBin.binIndexVec[i3+indexLength][i2+indexLength][i1+indexLength] = max_d_position;

				//this->allParams.oriBin.binIndexVec[i1+indexLength][i2+indexLength][i3+indexLength] = max_d_position;
				
				//this->allParams.oriBin.binIndexVec[i2+indexLength][i1+indexLength][i3+indexLength] = max_d_position;
			}
		}
	}

	int a1 = this->allParams.oriBin.binIndexVec[46][49][49];
	int a2 = this->allParams.oriBin.binIndexVec[49][46][49];
	
	std::vector<std::vector<double> > bin_centers(3, std::vector<double>(x.size(), 0));
	//bin_centers[0] = x;
	//bin_centers[1] = y;
	bin_centers[0] = y;
	bin_centers[1] = x;
	
	bin_centers[2] = z;
	this->allParams.oriBin.binCenters = bin_centers;
	
	std::vector<std::vector<double> > D(angleCount, std::vector<double>(angleCount, 0));
	std::vector<double> temp1(3, 0), temp2(3, 0);
	double acc = 0.0;
	for(int r = 0; r < angleCount; r++){
		for(int s = r+1; s < angleCount; s++){
			acc = 0.0;
			temp1[0] = bin_centers[0][r]; temp1[1] = bin_centers[1][r]; temp1[2] = bin_centers[2][r];
			temp2[0] = bin_centers[0][s]; temp2[1] = bin_centers[1][s]; temp2[2] = bin_centers[2][s];
			acc = std::inner_product(temp1.begin(), temp1.end(), temp2.begin(), acc); 

			D[r][s] = acc; D[s][r] = acc;
		}
	}
	
	// Using a map to sort numbers while preserving their index. The order for keys with same values might not be 
	// maintained after inserting in the map.
	int n_index_of_interest = this->allParams.oriBin.nLastIndicesOfInterest;
	VectorType2D nbr(angleCount, VectorType1D(n_index_of_interest, 0));

	/*VectorType1D D_sorted_index, D_sorted_index_of_interest;
	std::multimap<double, int> D_column_map;
	std::multimap<double, int>::iterator iter;
	for(int i = 0; i < D[0].size(); i++){
		
		for(int j = 0; j < D[0].size(); j++)
			D_column_map.insert(std::pair<double, int>(D[i][j], j));
		
		for(iter = D_column_map.begin(); iter != D_column_map.end(); iter++)
			D_sorted_index.push_back((*iter).second);
		
		std::reverse(D_sorted_index.begin(), D_sorted_index.end()); // for descending order
		
		nbr[i] = VectorType1D(D_sorted_index.begin(), D_sorted_index.begin() + n_index_of_interest);

		D_column_map.clear();
	}*/
	
	std::vector<arrayElement> D_col_vec;
	for(int i = 0; i < D.size(); i++){
	
		for(int j = 0; j < D[i].size(); j++)
			D_col_vec.push_back(arrayElement(D[i][j], j));
		
		std::stable_sort(D_col_vec.begin(), D_col_vec.end(), arrayElement::CompareArrayElementsDescendPreserveIndex);
		//std::stable_sort(D_col_vec.begin(), D_col_vec.end(), CompareArrayElementsPreserveIndex2);		
		//std::reverse(D_col_vec.begin(), D_col_vec.end());

		for(int j = 0; j < n_index_of_interest; j++)
			nbr[i][j] = D_col_vec[j].getIndex();
		
		D_col_vec.clear();
	}

	this->allParams.oriBin.nbr = nbr;	

	std::cout << "SphericalBinInfo computed." << std::endl;
}

void ftkVesselTracer::ComputeAllPrimaryNodes(void){
	
	this->globalStatsInput.volumeMax = Common::NormalizeData(this->inputData, this->normalizedInputData);
	this->ComputeSeeds();
	this->FitSphereAndSortNodes();
}

void ftkVesselTracer::ComputeSeeds(void){
	
	std::cout << "Computing initial seed points..." << std::endl;

	//ImageType3D::SizeType size = this->normalizedInputData->GetBufferedRegion().GetSize();
	ImageType3D::SizeType size = this->normalizedInputData->GetLargestPossibleRegion().GetSize();

	int grid_spacing = this->allParams.nodeDetectionParams.gridSpacing;

	//VolumeOfInterestFilterType::Pointer sub_volume_filter = VolumeOfInterestFilterType::New();
	//sub_volume_filter->SetInput(this->normalizedInputData);

	StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();

	ImageType3D::IndexType starting_index, max_index, min_index;
	ImageType3D::SizeType sub_volume_size; //, size1;
	ImageType3D::RegionType sub_volume_region;
	ImageType3D::Pointer sub_volume;
	//MinMaxCalculatorType::Pointer min_max_calculator = MinMaxCalculatorType::New();

	
	int d1 = size[0], d2 = size[1], d3 = size[2];
	//int j1 = 0, j2 = 0, j3 = 0;
	//double max_val = 0, min_val = 0, pixel1 = 0, pixel2 = 0, pixel3 = 0;
	int seed_count = 0;

	#pragma omp parallel for private(starting_index, max_index, min_index, sub_volume_size, sub_volume_region, sub_volume) 
	for(int i1 = 0; i1 < d1; i1 = i1 + grid_spacing){
		for(int i2 = 0; i2 < d2; i2 = i2 + grid_spacing){
			for(int i3 = 0; i3 < d3; i3 = i3 + grid_spacing){
			
				int j1 = std::min(i1 + grid_spacing - 1, d1);
				int j2 = std::min(i2 + grid_spacing - 1, d2);
				int j3 = std::min(i3 + grid_spacing - 1, d3);
				
				//starting_index[0] = i2; starting_index[1] = i1; starting_index[2] = i3;
				//sub_volume_size[0] = j2 - i2; sub_volume_size[1] = j1 - i1; sub_volume_size[2] = j3 - i3;
				starting_index[0] = i1; starting_index[1] = i2; starting_index[2] = i3;
				sub_volume_size[0] = j1 - i1; sub_volume_size[1] = j2 - i2; sub_volume_size[2] = j3 - i3;
				
				sub_volume_region.SetIndex(starting_index);
				sub_volume_region.SetSize(sub_volume_size);
				
				VolumeOfInterestFilterType::Pointer sub_volume_filter = VolumeOfInterestFilterType::New();
				sub_volume_filter->SetInput(this->normalizedInputData);
				sub_volume_filter->SetRegionOfInterest(sub_volume_region);
				sub_volume_filter->Update();
				sub_volume = sub_volume_filter->GetOutput();
				
				// Node statistics can be computed here to reject some seeds which lie in the background
				/*stats_filter->SetInput(sub_volume);
				stats_filter->Update();
				float mean = stats_filter->GetMean();
				float std = stats_filter->GetSigma();
				float max = stats_filter->GetMaximum();
				float min = stats_filter->GetMinimum();

				/*size1 = sub_volume->GetBufferedRegion().GetSize();
				double pixel1 = this->normalizedInputData->GetPixel(starting_index);
			    double pixel2 = sub_volume->GetPixel(starting_index);
				double pixel3 = this->inputData->GetPixel(starting_index);*/
				
				MinMaxCalculatorType::Pointer min_max_calculator = MinMaxCalculatorType::New();
				min_max_calculator->SetImage(sub_volume);
				min_max_calculator->Compute();
				
				double max_val = min_max_calculator->GetMaximum();
				double min_val = min_max_calculator->GetMinimum();
				max_index = min_max_calculator->GetIndexOfMaximum();
				min_index = min_max_calculator->GetIndexOfMinimum();

				//this->initialSeeds[seed_count].x = max_index[0];
				//this->initialSeeds[seed_count].y = max_index[1];
				//this->initialSeeds[seed_count].z = max_index[2];
				
				//max_index = min_index;	
				
				#pragma omp critical
				this->initialSeeds.push_back(Node(max_index[0] + i1, max_index[1] + i2, max_index[2] + i3, max_val));
				
				//#pragma omp critical
				//seed_count++;
			}
		}
	}
	
	std::cout << "Seed points computed." << std::endl;

	//visualizing the seed nodes
	//this->VisualizeNodesWithData3D(this->initialSeeds, false);
}

void ftkVesselTracer::VisualizeNodesWithData3D(std::vector<Node> node_vec, bool view_as_point){
	
	int n_seeds = node_vec.size();
	
	// Prepare the vtkVolume for rendering the data

	ITKToVTKConnectorType::Pointer ITK_to_VTK_connector = ITKToVTKConnectorType::New();

	ITK_to_VTK_connector->SetInput(this->inputDataForRendering);
	ITK_to_VTK_connector->Update();

	vtkSmartPointer<vtkImageData> vtk_image = ITK_to_VTK_connector->GetOutput();

	// Testing vtk image
	//vtk_image->PrintSelf(std::cout, vtkIndent(0));

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(1.0, 1.0, 1.0);
	
	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	renderer->SetActiveCamera(camera);

	vtkSmartPointer<vtkRenderWindow> render_window = vtkSmartPointer<vtkRenderWindow>::New();
	render_window->AddRenderer(renderer);
	
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	render_window_interactor->SetRenderWindow(render_window);
	
	vtkSmartPointer<vtkPiecewiseFunction> opacity_transfer_function = vtkSmartPointer<vtkPiecewiseFunction>::New();
	opacity_transfer_function->AddPoint(2, 0.0);
	opacity_transfer_function->AddPoint(10, 0.1);
	
	vtkSmartPointer<vtkColorTransferFunction> color_transfer_function = vtkSmartPointer<vtkColorTransferFunction>::New();
	color_transfer_function->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	color_transfer_function->AddRGBPoint((10 * 255)/100, 1, 1, 1); //blue
	color_transfer_function->AddRGBPoint((45 * 255)/100, 0, .01, 0); //green
	color_transfer_function->AddRGBPoint((150 * 255)/100, .01, 0, 0); //red
	
	vtkSmartPointer<vtkVolumeProperty> volume_property = vtkSmartPointer<vtkVolumeProperty>::New();
	volume_property->SetColor(color_transfer_function);
	volume_property->SetScalarOpacity(opacity_transfer_function);
	volume_property->ShadeOn();
	volume_property->SetInterpolationTypeToLinear();
	
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volume_mapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	
	volume_mapper->SetInput(vtk_image);
	volume_mapper->SetBlendModeToComposite();
	
	vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(volume_mapper); 
	volume->SetProperty(volume_property);
	volume->SetPosition(0, 0, 0);
	volume->SetPickable(0);
	volume->Update();

	// Prepare vtkActors for each node. Every node is displayed as a sphere

	std::vector<vtkSmartPointer<vtkSphereSource> > node_spheres;
	std::vector<vtkSmartPointer<vtkPolyDataMapper> > node_mappers;
	std::vector<vtkSmartPointer<vtkActor> > node_actors;

	for(int i = 0; i < n_seeds; i++){
		
		vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
		vtkSmartPointer<vtkPolyDataMapper> poly_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();

		sphere->SetCenter(node_vec[i].x, node_vec[i].y, node_vec[i].z);

		if(view_as_point == true)
			sphere->SetRadius(1);
		else
			sphere->SetRadius(node_vec[i].scale);
		
		poly_mapper->SetInputConnection(sphere->GetOutputPort());
		actor->SetMapper(poly_mapper);
		actor->GetProperty()->SetColor(0, 1, 0); //green
		
		node_spheres.push_back(sphere);
		node_mappers.push_back(poly_mapper);
		node_actors.push_back(actor);		
	}

	for(int i = 0; i < n_seeds; i++)		
		renderer->AddActor(node_actors[i]);

	renderer->AddVolume(volume);

	renderer->ResetCamera();
	
	render_window_interactor->Initialize();
	render_window->Render();
	render_window_interactor->Start();
}

Node::Node(){

	this->x = 0;
	this->y = 0;
	this->z = 0;
	this->intensity = 0;
	this->scale = 1;
	this->likelihood = 0;
	this->isValid = true;
	this->nHoodScale = this->scale;
	this->bandNhood = this->scale;
	this->bandArea = 0;
	this->foregroundArea = 0;
	this->backgroundArea = 0;
	this->meanBackgroundIntensity = 0;
	this->meanForegroundIntensity = 255;
	this->exitIter = 0;

	this->traceQuality = 0;
	this->parentIDLength = 0;
	this->isPrimary = false;
	this->isSecondary = false;
	this->nHoodSecondaryMultiplier = 0;
	this->nHoodScaleSecondary = 0;
	
	this->secondaryNodeSearchRad = 0;
	this->xInitSecondary = 0;
	this->yInitSecondary = 0;
	this->zInitSecondary = 0;

	
	this->ID = -1;
	this->NBranches = 0;
}

Node::Node(double x, double y, double z, PixelType intensity){

	this->x = x;
	this->y = y;
	this->z = z;
	this->intensity = intensity;

	this->scale = 1;
	this->likelihood = 0;
	this->isValid = true;
	this->nHoodScale = 2*this->scale;
	this->bandNhood = 0.5;
	this->bandArea = 0;
	this->foregroundArea = 0;
	this->backgroundArea = 0;
	this->meanBackgroundIntensity = 0;
	this->meanForegroundIntensity = 255;
	this->exitIter = 0;

	this->traceQuality = 0;
	this->parentIDLength = 4;
	this->isPrimary = true;
	this->isSecondary = false;
	this->nHoodSecondaryMultiplier = 2;
	this->nHoodScaleSecondary = this->nHoodSecondaryMultiplier * this->scale;
	
	this->secondaryNodeSearchRad = 0.0;
	this->xInitSecondary = 0;
	this->yInitSecondary = 0;
	this->zInitSecondary = 0;


	this->ID = -1;
	this->NBranches = 0;
}

bool compareNodesByLikelihood(Node n1, Node n2){	
	return (n1.likelihood < n2.likelihood); 
}

inline double Node::ComputeNorm(Node n){

	return(std::sqrt((n.x * n.x) + (n.y * n.y) + (n.z * n.z)));
}
void ftkVesselTracer::FitSphereAndSortNodes(void){

	std::cout << "Processing the seeds based on the model fits, likelihood and distance... " << std::endl;

	//seeds before fitting
	//this->VisualizeNodesWithData3D(this->initialSeeds, false);

	std::vector<Node> filteredPrimaryNodes;
	filteredPrimaryNodes.resize(this->initialSeeds.size());

	#pragma omp parallel for
	for(int i = 0; i < 50; i++){
	//for(int i = 0; i < this->initialSeeds.size(); i++){
		this->FitSphereAtNode(this->initialSeeds[i]);
		//std::cout << i << " Scale: " << this->initialSeeds[i].scale << " Likelihood: " << this->initialSeeds[i].likelihood << " Last iter: " << this->initialSeeds[i].exitIter << std::endl; 
		
		//#pragma omp critical
			if((this->initialSeeds[i].isValid == true) && (this->initialSeeds[i].likelihood > this->allParams.nodeDetectionParams.likelihoodThresholdPrimary)){
				//this->primaryNodes.push_back(this->initialSeeds[i]);
				filteredPrimaryNodes[i] = this->initialSeeds[i];

				//std::cout << " Scale: " << this->initialSeeds[i].scale << " Likelihood: " << this->initialSeeds[i].likelihood << std::endl;
			}
	}

	for(int i = 0; i < filteredPrimaryNodes.size(); i++){
		if(filteredPrimaryNodes[i].likelihood > 0)
			this->primaryNodes.push_back(filteredPrimaryNodes[i]);
	}

	//visualize the fitted models and the filtered nodes
	//this->VisualizeNodesWithData3D(this->initialSeeds, false);
	//this->VisualizeNodesWithData3D(this->primaryNodes, false);

	this->SortAndFilterPrimaryNodes();
	
	//final primary nodes
	this->VisualizeNodesWithData3D(this->primaryNodesAfterHitTest, false);

	std::cout << "Seed processing completed for primary nodes. " << std::endl;
	
	// Testing the model fitting for specific seeds
	//std::vector<Node> testing_seeds(this->initialSeeds.begin()+500, this->initialSeeds.begin()+510);
	//std::vector<Node> testing_seeds(this->initialSeeds.begin()+292, this->initialSeeds.begin()+293);
	/*std::vector<Node> testing_seeds(this->initialSeeds.begin()+200, this->initialSeeds.begin()+201);
	this->VisualizeNodesWithData3D(testing_seeds, false);
	for(int i = 0; i < testing_seeds.size(); i++)
		this->FitSphereAtNode(testing_seeds[i]);
	//std::cout << " Scale: " << testing_seeds[i].scale << " Likelihood: " << testing_seeds[i].likelihood << << " Position: " << seed.x << ", " << seed.y << ", " << seed.z " Last iter: " << testing_seeds[i].exitIter << std::endl; 
	this->VisualizeNodesWithData3D(testing_seeds, false);
	*/	
}

void ftkVesselTracer::FitSphereAtNode(Node& seed){

	//testing on a pseudo node
	// Same pixel apparently has different values in ITK and Matlab!!
	/*seed.x = 21;
	seed.y = 17;
	seed.z = 16;
	ImageType3D::IndexType test_index;
	test_index[0] = seed.x; test_index[1] = seed.y; test_index[2] = seed.z;
	PixelType test_pixel = this->inputData->GetPixel(test_index);
	PixelType test_pixel2 = this->normalizedInputData->GetPixel(test_index);
	*/
	
	// SET USING FUNCTION CALLS
	seed.scale = 4;
	seed.nHoodScale = 2 * seed.scale;
	seed.likelihood = 0;
	
	// Very important step!! 
	this->allParams.nodeDetectionParams.initByDefaultValues();

	int node_iter = this->allParams.nodeDetectionParams.iterNPrimaryNode;

	for(int i = 0; i < node_iter; i++){

		//std::cout << " Scale: " << seed.scale << " Position: " << seed.x << ", " << seed.y << ", " << seed.z << "  Likelihood: " << seed.likelihood << std::endl;
			
		this->UpdateAppearanceVectorized(seed);
		
		if(seed.isValid == false)
			break; 

		this->UpdateModel(seed, i);

		if(this->ExitModelFitting(seed, i) == true)
			break;
	}
}

void ftkVesselTracer::FitSphereAtNodeSecondary(Node& primary_node, Node& secondary_node, std::vector<double> dir_vec){

	Node anchor_node;
	anchor_node.x = primary_node.x; 
	anchor_node.y = primary_node.y; 
	anchor_node.z = primary_node.z;
	
	double anchor_dir_norm = Node::ComputeNorm(Node(dir_vec[0], dir_vec[1], dir_vec[2], 0.0));
	anchor_node.dirX.push_back(dir_vec[0]/anchor_dir_norm);
	anchor_node.dirY.push_back(dir_vec[1]/anchor_dir_norm);	
	anchor_node.dirZ.push_back(dir_vec[2]/anchor_dir_norm);

	anchor_node.secondaryNodeSearchRad = this->allParams.nodeDetectionParams.primaryNodeSearchRadFactor * primary_node.scale;
	
	anchor_node.xInitSecondary = primary_node.x + anchor_node.dirX[0] * anchor_node.secondaryNodeSearchRad;
	anchor_node.yInitSecondary = primary_node.y + anchor_node.dirY[0] * anchor_node.secondaryNodeSearchRad;
	anchor_node.zInitSecondary = primary_node.z + anchor_node.dirZ[0] * anchor_node.secondaryNodeSearchRad;


	secondary_node.isSecondary = true;
	secondary_node.x = anchor_node.xInitSecondary;
	secondary_node.y = anchor_node.yInitSecondary;
	secondary_node.z = anchor_node.zInitSecondary;
	secondary_node.scale = primary_node.scale;
	secondary_node.nHoodScale = 2.0 * secondary_node.scale;
	secondary_node.dirX.clear();
	secondary_node.dirY.clear();
	secondary_node.dirZ.clear();
	secondary_node.meanForegroundIntensity = primary_node.meanForegroundIntensity;
	secondary_node.meanBackgroundIntensity = primary_node.meanBackgroundIntensity;
	secondary_node.likelihood = 0.0;
	secondary_node.bandNhood = 0.5; // FUNCTION CALL REQUIRED
	secondary_node.sphHistRegionBased = primary_node.sphHistRegionBased;

	//Very important step !!
	this->allParams.nodeDetectionParams.initByDefaultValues();
	this->allParams.nodeDetectionParams.dtX = this->allParams.nodeDetectionParams.dtXSecondary;
	this->allParams.nodeDetectionParams.dtY = this->allParams.nodeDetectionParams.dtYSecondary;
	this->allParams.nodeDetectionParams.dtZ = this->allParams.nodeDetectionParams.dtZSecondary;
	this->allParams.nodeDetectionParams.dtScale = this->allParams.nodeDetectionParams.dtScaleSecondary;

	int node_iter = this->allParams.nodeDetectionParams.iterNPrimaryNode;

	for(int i = 0; i < node_iter; i++){

		//std::cout << " Scale: " << secondary_node.scale << " Position: " << secondary_node.x << ", " << secondary_node.y << ", " << secondary_node.z << "  Likelihood: " << secondary_node.likelihood << std::endl;
		
		this->UpdateAppearanceVectorized(secondary_node);

		if(secondary_node.isValid == false)
			break;
		
		this->UpdateModelSecondary(secondary_node, anchor_node, i);

		if(this->ExitModelFitting(secondary_node, i) == true)
			break;
	}
}

void ftkVesselTracer::UpdateAppearanceVectorized(Node& seed){
	
	ImageType3D::IndexType seed_index;
	seed_index[0] = seed.x; seed_index[1] = seed.y; seed_index[2] = seed.z;
	ImageType3D::SizeType size = this->normalizedInputData->GetLargestPossibleRegion().GetSize();

	//if(this->normalizedInputData->GetBufferedRegion().IsInside(seed_index) == false)
	if(this->normalizedInputData->GetLargestPossibleRegion().IsInside(seed_index) == false)
	//if((seed.x > 0) && (seed.x < size[0]) && (seed.y > 0) && (seed.y < size[1]) && (seed.z > 0) && (seed.z < size[2]))
		seed.isValid = false;
	
	ImageType3D::IndexType nhood_start_index, nhood_end_index;
	ImageType3D::SizeType nhood_size;
	ImageType3D::RegionType sub_vol;
	
	nhood_start_index[0] = floor(seed.x - seed.nHoodScale + 0.5); 
	nhood_start_index[1] = floor(seed.y - seed.nHoodScale + 0.5); 
	nhood_start_index[2] = floor(seed.z - seed.nHoodScale + 0.5);
	nhood_end_index[0] = floor(seed.x + seed.nHoodScale + 0.5);
	nhood_end_index[1] = floor(seed.y + seed.nHoodScale + 0.5);
	nhood_end_index[2] = floor(seed.z + seed.nHoodScale + 0.5);

	nhood_size[0] = 2*seed.nHoodScale; nhood_size[1] = 2*seed.nHoodScale; nhood_size[2] = 2*seed.nHoodScale;
	
	/*VolumeOfInterestFilterType::Pointer sub_vol_filter = VolumeOfInterestFilterType::New();
	sub_vol_filter->SetInput(this->normalizedInputData);
	sub_vol.SetIndex(nhood_start_index);
	sub_vol.SetSize(nhood_size);
	sub_vol_filter->SetRegionOfInterest(sub_vol);
	sub_vol_filter->Update();
	*/
	
	seed.bandArea = 0;
	seed.foregroundArea = 0;
	seed.backgroundArea = 0;
	seed.meanForegroundIntensity = 0;
	seed.meanBackgroundIntensity = 0;
	seed.xNormalizedInBand.clear();
	seed.yNormalizedInBand.clear();
	seed.zNormalizedInBand.clear();
	seed.intensityInBand.clear();
	seed.gxInBand.clear();
	seed.gyInBand.clear();
	seed.gzInBand.clear();

	ImageType3D::IndexType current_index;
	PixelType x_normalized = 0, y_normalized = 0, z_normalized = 0, radial_dist = 0;
	bool is_inside_volume = true, is_outside_volume = false, is_on_foreground = true, is_on_background = false, contributes_to_energy = true, energy_band_outside_volume = false;
	int in_volume_count = 0, out_volume_count = 0;
	std::vector<double> foreground_values, background_values;
	
	for(int k = nhood_start_index[2]; k <= nhood_end_index[2]; k++){
		for(int i = nhood_start_index[0]; i <= nhood_end_index[0]; i++){
			for(int j = nhood_start_index[1]; j <= nhood_end_index[1]; j++){
			//for(int k = nhood_start_index[2]; k <= nhood_end_index[2]; k++){
				
				current_index[0] = i; current_index[1] = j; current_index[2] = k; 
				
				//if(this->normalizedInputData->GetBufferedRegion().IsInside(current_index) == true){
				if(this->normalizedInputData->GetLargestPossibleRegion().IsInside(current_index) == true){
				//if((i > 0) && (i < size[0]) && (j > 0) && (j < size[1]) && (k > 0) && (k < size[2])){
					is_inside_volume = true;
					in_volume_count++;
				}
				else{
					is_inside_volume = false;
					out_volume_count++;
				}
								
				x_normalized = (i - seed.x)/seed.scale;
				y_normalized = (j - seed.y)/seed.scale;
				z_normalized = (k - seed.z)/seed.scale;
				radial_dist = std::sqrt((x_normalized*x_normalized) + (y_normalized*y_normalized) + (z_normalized*z_normalized));

				if((radial_dist <= 1.0) && (is_inside_volume == true))
				//if(((radial_dist < 1.0) || (std::abs(radial_dist - 1.0) < 0.000001)) && (is_inside_volume == true))
					is_on_foreground = true;
				else
					is_on_foreground = false;

				if((is_on_foreground == false) && (is_inside_volume == true))
					is_on_background = true;
				else
					is_on_background = false;

				if((abs(radial_dist - 1) < seed.bandNhood) && (is_inside_volume == true))		
					contributes_to_energy = true;
				else
					contributes_to_energy = false;

				if((abs(radial_dist - 1) < seed.bandNhood) && (is_inside_volume == false))		
					energy_band_outside_volume = true;
				else
					energy_band_outside_volume = false;

				if(contributes_to_energy == true){
					seed.xNormalizedInBand.push_back(x_normalized/radial_dist);
					seed.yNormalizedInBand.push_back(y_normalized/radial_dist);
					seed.zNormalizedInBand.push_back(z_normalized/radial_dist);
					seed.intensityInBand.push_back(this->normalizedInputData->GetPixel(current_index));
					seed.gxInBand.push_back(this->gx->GetPixel(current_index));
					seed.gyInBand.push_back(this->gy->GetPixel(current_index));
					seed.gzInBand.push_back(this->gz->GetPixel(current_index));

					seed.bandArea++;
				}

				if(energy_band_outside_volume == true){
					seed.xNormalizedInBand.push_back(x_normalized/radial_dist);
					seed.yNormalizedInBand.push_back(y_normalized/radial_dist);
					seed.zNormalizedInBand.push_back(z_normalized/radial_dist);
					seed.intensityInBand.push_back(0);
					seed.gxInBand.push_back(0);
					seed.gyInBand.push_back(0);
					seed.gzInBand.push_back(0);
					
					seed.bandArea++;
				}
				
				if(is_on_foreground == true){

					//std::cout << "Position: " << i << ", " << j << ", " << k << " Value: " << this->normalizedInputData->GetPixel(current_index) << std::endl;
					
					foreground_values.push_back(this->normalizedInputData->GetPixel(current_index));
					
					seed.meanForegroundIntensity = seed.meanForegroundIntensity + this->normalizedInputData->GetPixel(current_index);
					seed.foregroundArea++;
				}
				if(is_on_background == true){

					background_values.push_back(this->normalizedInputData->GetPixel(current_index));

					seed.meanBackgroundIntensity = seed.meanBackgroundIntensity + this->normalizedInputData->GetPixel(current_index);
					seed.backgroundArea++;
				}

			} 
		}
	}

	if(seed.isSecondary == true){
		
		if(foreground_values.size() != 0){
			// Calculate the median
			std::sort(foreground_values.begin(), foreground_values.end());

			seed.meanForegroundIntensity = *(foreground_values.begin() + foreground_values.size() / 2);
		}
		else
			seed.meanForegroundIntensity = 0.0;

		if(background_values.size() != 0){
			std::sort(background_values.begin(), background_values.end());
			seed.meanBackgroundIntensity = *(background_values.begin() + background_values.size() / 2);
		}
		else
			seed.meanBackgroundIntensity = 0.0;
	}
	else{

		// Calculate the mean
		seed.meanForegroundIntensity = seed.meanForegroundIntensity / (double)seed.foregroundArea;
		seed.meanBackgroundIntensity = seed.meanBackgroundIntensity / (double)seed.backgroundArea;
	}

	seed.likelihood = seed.meanForegroundIntensity - seed.meanBackgroundIntensity;

	// If likelood is less than 0.001, do what? According to the thesis, the subvolume should be discarded. But here, the likelihood value is
	// increased if it is less then 0.001
	/*double mean_of_means = 0.0;
	if(seed.isSecondary == true && seed.likelihood < this->allParams.nodeDetectionParams.increaseLikelihoodThreshold){
		
		mean_of_means = (seed.meanForegroundIntensity + seed.meanBackgroundIntensity) / 2.0;
		seed.meanForegroundIntensity = mean_of_means + this->allParams.nodeDetectionParams.discardNodeLikelihoodThreshold;
		seed.meanBackgroundIntensity = mean_of_means - this->allParams.nodeDetectionParams.discardNodeLikelihoodThreshold;
		if(seed.meanBackgroundIntensity < 0)
			seed.meanBackgroundIntensity = 0;
		seed.likelihood = seed.meanForegroundIntensity - seed.meanBackgroundIntensity;
	}*/
	
	// Set seed as invalid if majority of its neighbourhood lies outside the data limits
	if(in_volume_count < out_volume_count)
		seed.isValid = false;
	else
		seed.isValid = true;	
}

void ftkVesselTracer::UpdateModel(Node& seed, int iter_number){

	double last_x = seed.x, last_y = seed.y, last_z = seed.z;
	double last_scale = seed.scale;

	std::vector<double> inside_region_term, outside_region_term;
	for(int i = 0; i < seed.intensityInBand.size(); i++){
		inside_region_term.push_back(std::abs(seed.meanForegroundIntensity - seed.intensityInBand[i]));
		outside_region_term.push_back(std::abs(PixelType(seed.meanBackgroundIntensity - seed.intensityInBand[i])));
	}
	
	std::vector<double> del_energy, region_based_term, gvf_based_term;
	//std::transform(inside_region_term.begin(), inside_region_term.end(), outside_region_term.begin(), del_energy.begin(), std::minus<double>());
	for(int i = 0; i < inside_region_term.size(); i++){
		region_based_term.push_back(inside_region_term[i] - outside_region_term[i]);
	}
	
	double region_weight = this->allParams.nodeDetectionParams.regionBasedTermWeight;
	double gvf_weight = this->allParams.nodeDetectionParams.edgeBasedTermWeight;
	if(iter_number < this->allParams.nodeDetectionParams.iterNForOnlyRegionBasedTerm)
		del_energy = region_based_term;
	else{
		for(int i = 0; i < seed.bandArea; i++)
			gvf_based_term.push_back(seed.xNormalizedInBand[i] * seed.gxInBand[i] + seed.yNormalizedInBand[i] * seed.gyInBand[i] + seed.zNormalizedInBand[i] * seed.gzInBand[i]);
		
		for(int i = 0; i < seed.bandArea; i++)
			del_energy.push_back(region_weight * region_based_term[i] - gvf_weight * gvf_based_term[i]);
	}
	
	double acc_energy = std::accumulate(del_energy.begin(), del_energy.end(), 0.0);
	double dscale = acc_energy / seed.bandArea;

	//std::cout << " Energy: " << acc_energy << std::endl;
	
	double dx = 0, dy = 0, dz = 0;
	for(int i = 0; i < seed.bandArea; i++){
		dx = dx + (seed.xNormalizedInBand[i] * del_energy[i]/seed.bandArea);
		dy = dy + (seed.yNormalizedInBand[i] * del_energy[i]/seed.bandArea);
		dz = dz + (seed.zNormalizedInBand[i] * del_energy[i]/seed.bandArea);
	}
	
	dscale = dscale * this->allParams.nodeDetectionParams.dtScale;
	dx = dx * this->allParams.nodeDetectionParams.dtX;
	dy = dy * this->allParams.nodeDetectionParams.dtY;
	dz = dz * this->allParams.nodeDetectionParams.dtZ;

	dscale = std::min(1.0, std::max(-1.0, dscale));
	dx = std::min(1.0, std::max(-1.0, dx));
	dy = std::min(1.0, std::max(-1.0, dy));
	dz = std::min(1.0, std::max(-1.0, dz));

	seed.scale = seed.scale - dscale;
	seed.x = seed.x - dx;
	seed.y = seed.y - dy;
	seed.z = seed.z - dz;

	if(seed.scale > this->allParams.nodeDetectionParams.maxVesselWidth)
		seed.scale = this->allParams.nodeDetectionParams.maxVesselWidth;
	if(seed.scale < this->allParams.nodeDetectionParams.minVesselWidth)
		seed.scale = this->allParams.nodeDetectionParams.minVesselWidth;

	if(this->GetSign(dx) != this->GetSign(this->allParams.nodeDetectionParams.dirX))
		this->allParams.nodeDetectionParams.dtX = this->allParams.nodeDetectionParams.dtX * this->allParams.nodeDetectionParams.primaryReversePositionRate;
	if(this->GetSign(dy) != this->GetSign(this->allParams.nodeDetectionParams.dirY))
		this->allParams.nodeDetectionParams.dtY = this->allParams.nodeDetectionParams.dtY * this->allParams.nodeDetectionParams.primaryReversePositionRate;
	if(this->GetSign(dz) != this->GetSign(this->allParams.nodeDetectionParams.dirZ))
		this->allParams.nodeDetectionParams.dtZ = this->allParams.nodeDetectionParams.dtZ * this->allParams.nodeDetectionParams.primaryReversePositionRate;

	this->allParams.nodeDetectionParams.dirX = this->GetSign(dx);
	this->allParams.nodeDetectionParams.dirY = this->GetSign(dy);
	this->allParams.nodeDetectionParams.dirZ = this->GetSign(dz);

	this->allParams.nodeDetectionParams.chX[this->allParams.nodeDetectionParams.currentMonitoredIter] = seed.x - last_x;
	this->allParams.nodeDetectionParams.chY[this->allParams.nodeDetectionParams.currentMonitoredIter] = seed.y - last_y;
	this->allParams.nodeDetectionParams.chZ[this->allParams.nodeDetectionParams.currentMonitoredIter] = seed.z - last_z;

	
	if((dscale * this->allParams.nodeDetectionParams.dirScale) <= 0)
		this->allParams.nodeDetectionParams.dtScale = this->allParams.nodeDetectionParams.dtScale * this->allParams.nodeDetectionParams.primaryReverseScaleRate;

	this->allParams.nodeDetectionParams.dirScale = this->GetSign(dscale);
	this->allParams.nodeDetectionParams.chScale[this->allParams.nodeDetectionParams.currentMonitoredIter] = seed.scale - last_scale;
}

void ftkVesselTracer::UpdateModelSecondary(Node& seed, Node& anchor_node, int iter_number){

	double last_x = seed.x, last_y = seed.y, last_z = seed.z;
	double last_scale = seed.scale;
	
	std::vector<double> inside_region_term, outside_region_term;
	for(int i = 0; i < seed.intensityInBand.size(); i++){
		inside_region_term.push_back(std::abs(seed.meanForegroundIntensity - seed.intensityInBand[i]));
		outside_region_term.push_back(std::abs(PixelType(seed.meanBackgroundIntensity - seed.intensityInBand[i])));
	}
	
	std::vector<double> del_energy, region_based_term, gvf_based_term;
	//std::transform(inside_region_term.begin(), inside_region_term.end(), outside_region_term.begin(), del_energy.begin(), std::minus<double>());
	for(int i = 0; i < inside_region_term.size(); i++){
		region_based_term.push_back(inside_region_term[i] - outside_region_term[i]);
	}
	
	double region_weight = this->allParams.nodeDetectionParams.regionBasedTermWeightSecondary;
	double gvf_weight = this->allParams.nodeDetectionParams.edgeBasedTermWeightSecondary;
	if(iter_number < this->allParams.nodeDetectionParams.iterNForOnlyRegionBasedTerm)
		del_energy = region_based_term;
	else{
		for(int i = 0; i < seed.bandArea; i++)
			gvf_based_term.push_back(seed.xNormalizedInBand[i] * seed.gxInBand[i] + seed.yNormalizedInBand[i] * seed.gyInBand[i] + seed.zNormalizedInBand[i] * seed.gzInBand[i]);
		
		for(int i = 0; i < seed.bandArea; i++)
			del_energy.push_back(region_weight * region_based_term[i] - gvf_weight * gvf_based_term[i]);
	}

	double acc_energy = std::accumulate(del_energy.begin(), del_energy.end(), 0.0);
	double dscale = acc_energy / seed.bandArea;

	//std::cout << " Energy: " << acc_energy << std::endl;
	
	double dx = 0, dy = 0, dz = 0, dc = 0;
	for(int i = 0; i < seed.bandArea; i++){
		dx = dx + (seed.xNormalizedInBand[i] * del_energy[i]/seed.bandArea);
		dy = dy + (seed.yNormalizedInBand[i] * del_energy[i]/seed.bandArea);
		dz = dz + (seed.zNormalizedInBand[i] * del_energy[i]/seed.bandArea);
	}
	
	dscale = dscale * this->allParams.nodeDetectionParams.dtScale;
	dx = dx * this->allParams.nodeDetectionParams.dtX;
	dy = dy * this->allParams.nodeDetectionParams.dtY;
	dz = dz * this->allParams.nodeDetectionParams.dtZ;

	dscale = std::min(0.5, std::max(-0.5, dscale));
	dx = std::min(0.5, std::max(-0.5, dx));
	dy = std::min(0.5, std::max(-0.5, dy));
	dz = std::min(0.5, std::max(-0.5, dz));

	dc = anchor_node.dirX[0]*dx + anchor_node.dirY[0]*dy + anchor_node.dirZ[0]*dz;

	seed.scale = seed.scale - dscale;
	seed.x = seed.x - dx + dc*anchor_node.dirX[0];
	seed.y = seed.y - dy + dc*anchor_node.dirY[0];
	seed.z = seed.z - dz + dc*anchor_node.dirZ[0];

	if(seed.scale > this->allParams.nodeDetectionParams.maxVesselWidth)
		seed.scale = this->allParams.nodeDetectionParams.maxVesselWidth;
	if(seed.scale < this->allParams.nodeDetectionParams.minVesselWidth)
		seed.scale = this->allParams.nodeDetectionParams.minVesselWidth;

	std::vector<double> ch_dir(3, 0), ch_dir_normalized(3, 0);
	ch_dir[0] = seed.x - anchor_node.xInitSecondary;
	ch_dir[1] = seed.y - anchor_node.yInitSecondary;
	ch_dir[2] = seed.z - anchor_node.zInitSecondary;
	double ch_norm = Node::ComputeNorm(Node(ch_dir[0], ch_dir[1], ch_dir[2], 0));
	ch_dir_normalized[0] = ch_dir[0]/ch_norm;
	ch_dir_normalized[1] = ch_dir[1]/ch_norm;
	ch_dir_normalized[2] = ch_dir[2]/ch_norm;
	
	// Restrict the search space for secondary nodes
	if(ch_norm > this->allParams.nodeDetectionParams.secondarySearchConstraint * seed.scale){	
		seed.x = anchor_node.xInitSecondary + ch_dir[0] * (this->allParams.nodeDetectionParams.secondarySearchConstraint * seed.scale / ch_norm);
		seed.y = anchor_node.yInitSecondary + ch_dir[1] * (this->allParams.nodeDetectionParams.secondarySearchConstraint * seed.scale / ch_norm);
		seed.z = anchor_node.zInitSecondary + ch_dir[2] * (this->allParams.nodeDetectionParams.secondarySearchConstraint * seed.scale / ch_norm);
	}

	if(this->GetSign(dx) != this->GetSign(this->allParams.nodeDetectionParams.dirX))
		this->allParams.nodeDetectionParams.dtX = this->allParams.nodeDetectionParams.dtX * this->allParams.nodeDetectionParams.secondaryReversePositionRate;
	if(this->GetSign(dy) != this->GetSign(this->allParams.nodeDetectionParams.dirY))
		this->allParams.nodeDetectionParams.dtY = this->allParams.nodeDetectionParams.dtY * this->allParams.nodeDetectionParams.secondaryReversePositionRate;
	if(this->GetSign(dz) != this->GetSign(this->allParams.nodeDetectionParams.dirZ))
		this->allParams.nodeDetectionParams.dtZ = this->allParams.nodeDetectionParams.dtZ * this->allParams.nodeDetectionParams.secondaryReversePositionRate;

	this->allParams.nodeDetectionParams.dirX = this->GetSign(dx);
	this->allParams.nodeDetectionParams.dirY = this->GetSign(dy);
	this->allParams.nodeDetectionParams.dirZ = this->GetSign(dz);

	this->allParams.nodeDetectionParams.chX[this->allParams.nodeDetectionParams.currentMonitoredIter] = seed.x - last_x;
	this->allParams.nodeDetectionParams.chY[this->allParams.nodeDetectionParams.currentMonitoredIter] = seed.y - last_y;
	this->allParams.nodeDetectionParams.chZ[this->allParams.nodeDetectionParams.currentMonitoredIter] = seed.z - last_z;

	
	if((dscale * this->allParams.nodeDetectionParams.dirScale) <= 0)
		this->allParams.nodeDetectionParams.dtScale = this->allParams.nodeDetectionParams.dtScale * this->allParams.nodeDetectionParams.secondaryReverseScaleRate;

	this->allParams.nodeDetectionParams.dirScale = this->GetSign(dscale); // = dscale; // As given in Matlab code
	this->allParams.nodeDetectionParams.chScale[this->allParams.nodeDetectionParams.currentMonitoredIter] = seed.scale - last_scale;
}

int inline ftkVesselTracer::GetSign(double value){
	if(value < 0)
		return -1;
	if(value > 0) 
		return 1;
	return 0;
}

bool ftkVesselTracer::ExitModelFitting(Node& seed, int iter_number){
	
	bool exit_fitting = false;
	this->allParams.nodeDetectionParams.currentMonitoredIter = (iter_number % this->allParams.nodeDetectionParams.iterNMonitorParamChange); // + 1;
	
	if(seed.scale >= this->allParams.nodeDetectionParams.maxVesselWidth){
		seed.exitIter = iter_number;
		exit_fitting = true;
	}
	
	double total_ch_scale = 0, total_ch_x = 0, total_ch_y = 0, total_ch_z = 0, total_ch_position = 0;
	int iter_minimum = 0;
	if(seed.isSecondary == true)
		iter_minimum = this->allParams.nodeDetectionParams.iterNMinimumSecondary;
	else
		iter_minimum = this->allParams.nodeDetectionParams.iterNMinimum;

	if(iter_number > iter_minimum){
		for(int i = 0; i < this->allParams.nodeDetectionParams.iterNMonitorParamChange; i++){
			total_ch_scale = total_ch_scale + std::abs(this->allParams.nodeDetectionParams.chScale[i]);
			total_ch_x = total_ch_x + std::abs(this->allParams.nodeDetectionParams.chX[i]);
			total_ch_y = total_ch_y + std::abs(this->allParams.nodeDetectionParams.chY[i]);
			total_ch_z = total_ch_z + std::abs(this->allParams.nodeDetectionParams.chZ[i]);
		}
		total_ch_position = total_ch_x + total_ch_y + total_ch_z;

		if(total_ch_scale < this->allParams.nodeDetectionParams.minimumAccumulatedParamChange && total_ch_position < this->allParams.nodeDetectionParams.minimumAccumulatedParamChange){
			seed.exitIter = iter_number;
			exit_fitting = true;
		}
	}
	if(iter_number == this->allParams.nodeDetectionParams.iterNPrimaryNode)
		seed.exitIter = iter_number;
	
	return exit_fitting;
}

void ftkVesselTracer::SortAndFilterPrimaryNodes(void){

	// Remove nodes lying outside the image space
	ImageType3D::IndexType index;
	for(int i = 0; i < this->primaryNodes.size(); i++){
		
		index[0] = this->primaryNodes[i].x; index[1] = this->primaryNodes[i].y; index[2] = this->primaryNodes[i].z;

		if(this->normalizedInputData->GetLargestPossibleRegion().IsInside(index) == false)
			this->primaryNodes.erase(this->primaryNodes.begin() + i);
	}

	// Sort the primary nodes in descending order of likelihood
	std::sort(this->primaryNodes.begin(), this->primaryNodes.end(), compareNodesByLikelihood);
	std::reverse(this->primaryNodes.begin(), this->primaryNodes.end());

	// Remove nodes which are very close to each other (node with higher likelihood remains)
	this->primaryNodesAfterHitTest.push_back(this->primaryNodes[0]);
	double norm = 0;
	bool hit = false;
	
	for(int i = 1; i < this->primaryNodes.size(); i++){		
		
		Node current_primary_node = this->primaryNodes[i];
		hit = false;
		for(int j = 0; j < this->primaryNodesAfterHitTest.size(); j++){
			
			Node current_primary_node_final = this->primaryNodesAfterHitTest[j];
			norm = Node::ComputeNorm(Node(current_primary_node.x - current_primary_node_final.x, current_primary_node.y - current_primary_node_final.y, 
				current_primary_node.z - current_primary_node_final.z, 0));
			if(norm < (this->allParams.nodeDetectionParams.distanceThresholdPrimary * (current_primary_node.scale + current_primary_node_final.scale))){
				hit = true;
				break;
			}
		}
		if(hit == false)
			this->primaryNodesAfterHitTest.push_back(current_primary_node);
	}

	std::cout << "Final primary nodes: " << this->primaryNodesAfterHitTest.size() << std::endl;
	for(int i = 0; i < this->primaryNodesAfterHitTest.size(); i++)
		std::cout << " Scale: " << this->primaryNodesAfterHitTest[i].scale << " Likelihood: " << this->primaryNodesAfterHitTest[i].likelihood << std::endl;
}

void ftkVesselTracer::ComputeAllSecondaryNodes(void){

	// For testing purposes
	ImageType3D::IndexType test_index, test_index1;
	test_index[0] = 128; test_index[1] = 7; test_index[2] = 25;
	test_index1[0] = 6; test_index1[1] = 127; test_index1[2] = 24;
	PixelType test_pixel = this->inputData->GetPixel(test_index);
	PixelType test_pixel1 = this->inputData->GetPixel(test_index1);

	//this->primaryNodesAfterHitTest.erase(this->primaryNodesAfterHitTest.begin()+1, this->primaryNodesAfterHitTest.end());

	ImageType3D::IndexType test_index2;
	test_index2[0] = this->primaryNodesAfterHitTest[0].x; 
	test_index2[1] = this->primaryNodesAfterHitTest[0].y;
	test_index2[2] = this->primaryNodesAfterHitTest[0].z;
	PixelType test_pixel2 = this->inputData->GetPixel(test_index2);
	PixelType test_pixel3 = this->normalizedInputData->GetPixel(test_index2);
	
	
	clock_t start_tracing_time, stop_tracing_time;
	double tracing_time = 0;

	assert((start_tracing_time = clock()) != -1);
	
	/*Node a_node;
	a_node.scale = 3.8346;
	a_node.y = 124.7428;
	a_node.x = 64.4287;
	a_node.z = 19.1630;
	a_node.likelihood = 0.4169;
	a_node.meanForegroundIntensity = 0.5671;
	a_node.meanBackgroundIntensity = 0.1502;
	a_node.parentIDLength = 4;
	this->allParams.nodeDetectionParams.traceLengthCost = 0.5;
	a_node.nHoodSecondaryMultiplier = 2;
	this->primaryNodesAfterHitTest.erase(this->primaryNodesAfterHitTest.begin()+1, this->primaryNodesAfterHitTest.end());
	this->primaryNodesAfterHitTest[0] = a_node; */
	this->primaryNodesAfterHitTest.erase(this->primaryNodesAfterHitTest.begin()+1, this->primaryNodesAfterHitTest.end());
	
	//this->VisualizeNodesWithData3D(this->primaryNodesAfterHitTest, false);

	PriorityQueueType node_queue(compareNodes(2));
	std::vector<queue_element> light_node_queue(this->allParams.nodeDetectionParams.maxQueueSize, queue_element(0, this->allParams.nodeDetectionParams.infTraceQuality));
	std::vector<double> quality_array;	

	//Node a_node;
	double trace_quality = 0.0;
	for(int i = 0; i < this->primaryNodesAfterHitTest.size(); i++){
		//Node a_node = this->primaryNodesAfterHitTest[i];           

		// Quality estimate for single nodes. Nodes with high likelihood get low quality, nodes with low likelihood get high quality
		if(this->primaryNodesAfterHitTest[i].likelihood <= 0){
			this->primaryNodesAfterHitTest[i].traceQuality = this->allParams.nodeDetectionParams.maxTraceCost; //.traceQualityThreshold;
			quality_array.push_back(this->allParams.nodeDetectionParams.maxTraceCost);
			light_node_queue[i] = queue_element(i, this->allParams.nodeDetectionParams.maxTraceCost);
		}
		else{
			trace_quality = -1 * std::log(this->primaryNodesAfterHitTest[i].likelihood) + this->allParams.nodeDetectionParams.traceLengthCost;
			this->primaryNodesAfterHitTest[i].traceQuality = trace_quality;
			quality_array.push_back(trace_quality);
			light_node_queue[i] = queue_element(i, trace_quality);
		}
	
		this->primaryNodesAfterHitTest[i].parentID = std::vector<double>(this->primaryNodesAfterHitTest[i].parentIDLength, -1);
		//node_queue.push(this->primaryNodesAfterHitTest[i]);
	}

	//std::cout << "Primary node trace quality.. " << std::endl;
	//for(int i = 0; i < quality_array.size(); i++)
	//	std::cout << quality_array[i] << std::endl;
	
	int total_nodes_counter = -1, hit = 0, primary_counter = 0, hit_counter = 0, queue_iter = -1;
	int queue_size = this->primaryNodesAfterHitTest.size();
	Node current_node, dir_node;
	double dirX = 0, dirY = 0, dirZ = 0, norm = 0;
	std::vector<double> dir_hist; 
	std::vector<Node> badNodes, veryBadNodes;
	queue_element current_queue_element(0, 0);
	
	//while(queue_iter <= this->primaryNodesAfterHitTest.size() && queue_iter < this->allParams.nodeDetectionParams.maxQueueIter && !node_queue.empty()){
	//while(!node_queue.empty() && queue_iter < this->allParams.nodeDetectionParams.maxQueueIter){
	while(queue_iter <= queue_size && queue_iter < this->allParams.nodeDetectionParams.maxQueueIter){
		
		//current_node = node_queue.top();
		//node_queue.pop();
		
		this->GetBestTrace(light_node_queue, current_queue_element);

		// Break if queue is empty
		if(current_queue_element.second == this->allParams.nodeDetectionParams.infTraceQuality)
			break;

		current_node = this->primaryNodesAfterHitTest[current_queue_element.first];

		if(current_node.parentID[0] == -1)
			current_node.isPrimary = true;
		else
			current_node.isPrimary = false;
	
		// If the quality of trace is above the threshold, these nodes have low quality (ignore them)
		if(current_queue_element.second > this->allParams.nodeDetectionParams.traceQualityThreshold && current_node.isPrimary == false)
			continue;
		
		hit = this->TraceHitTest(current_node);

		if(hit == 2 || (hit == 1 && current_node.isPrimary == true))
			continue;

		if(current_node.isPrimary == true && hit == 0)
			primary_counter++;
		
		total_nodes_counter++;
		this->allNodes.push_back(current_node);

		if(hit == 1){
			hit_counter++;
			continue;
		} 
		
		if(current_node.isPrimary == true){
			
			dir_hist = std::vector<double>(this->allParams.oriBin.angleCount, 1.0/this->allParams.oriBin.angleCount);
			current_node.dirX.push_back(0);
			current_node.dirY.push_back(0);
			current_node.dirZ.push_back(0);
		}
		else{		
			dirX = this->allNodes[current_node.parentID[0]].x - current_node.x;
			dirY = this->allNodes[current_node.parentID[0]].y - current_node.y;
			dirZ = this->allNodes[current_node.parentID[0]].z - current_node.z;
			
			dir_node = Node(dirX, dirY, dirZ, 0);
			
			norm = Node::ComputeNorm(dir_node);
			dir_node = Node(dirX/norm, dirY/norm, dirZ/norm, 0);
			dir_hist = current_node.sphHistRegionBased;

			this->allNodes[total_nodes_counter].sphHistRegionBased.clear();
		}

		this->ComputeSecondaryNodeDirections(current_node, dir_hist);

		if(current_node.dirX.empty() || current_node.dirY.empty() || current_node.dirZ.empty())
			continue;

		std::vector<double> current_dir(3, 0);
		//Node secondary_node;
		for(int i = 0; i < current_node.dirX.size(); i++){
			
			// Only branches less than maxBranchAngle are considered
			if((current_node.dirX[i] * dirX + current_node.dirY[i] * dirY + current_node.dirZ[i] * dirZ) > std::cos(this->allParams.nodeDetectionParams.maxBranchAngle * vnl_math::pi / 180))
				continue;
			
			current_dir[0] = current_node.dirX[i];
			current_dir[1] = current_node.dirY[i];
			current_dir[2] = current_node.dirZ[i];
			
			Node secondary_node;
			secondary_node.parentIDLength = current_node.parentIDLength;
			secondary_node.parentID = std::vector<double>(current_node.parentIDLength, -1);

			this->FitSphereAtNodeSecondary(current_node, secondary_node, current_dir);			
			
			if(secondary_node.isValid == false){ // || secondary_node.likelihood <= 0){
				badNodes.push_back(secondary_node);
				quality_array.push_back(this->allParams.nodeDetectionParams.maxTraceCost);
				continue;
			}
			if(secondary_node.likelihood <= 0.0){
				veryBadNodes.push_back(secondary_node);
				quality_array.push_back(2*this->allParams.nodeDetectionParams.maxTraceCost);
			}
			
			// Parent IDs of secondary node = circshift(parent IDs of current_node)
			secondary_node.parentID[0] = current_node.parentID.back(); //*(current_node.parentID.begin() + current_node.parentID.size()-1);
			for(int j = 1; j < current_node.parentID.size(); j++){
				secondary_node.parentID[j] = current_node.parentID[j-1];
			}
			secondary_node.parentID[0] = total_nodes_counter;
			
			//queue_size++;
			this->primaryNodesAfterHitTest.push_back(secondary_node);

			trace_quality = this->computeTraceQuality(secondary_node);
	
			std::cout << "Trace quality: " << trace_quality << std::endl;
			
			//node_queue.push(current_node);
			//node_queue.push(secondary_node);

			quality_array.push_back(trace_quality);

			light_node_queue[queue_size] = queue_element(queue_size, trace_quality);
			queue_size++;
		}
		queue_iter++;
		
		//this->VisualizeNodesWithData3D(this->primaryNodesAfterHitTest, false);
		//this->VisualizeNodesWithData3D(this->allNodes, false);
	}

	stop_tracing_time = clock();
	tracing_time = (double)(stop_tracing_time - start_tracing_time)/CLOCKS_PER_SEC;

	std::cout << "Tracing took " << tracing_time << " seconds." << std::endl;

	//Visualize all nodes
	//this->VisualizeNodesWithData3D(this->allNodes, true);
	this->VisualizeNodesWithData3D(this->allNodes, false);

	std::vector<Node> nodesWithNonPositiveLikelihood;
	for(int i = 0; i < this->allNodes.size(); i++){
		if(this->allNodes[i].likelihood <= 0.0)
			nodesWithNonPositiveLikelihood.push_back(this->allNodes[i]);
	}
	this->VisualizeNodesWithData3D(nodesWithNonPositiveLikelihood, false);

	if(badNodes.empty() == false)
		this->VisualizeNodesWithData3D(badNodes, false);

	//Write all nodes to a file
	this->writeNodesToFile(this->allNodes, std::string("AllNodes.txt"));

	//write quality array
	ofstream nodes_file_stream;
	nodes_file_stream.open("QualityArray.txt", std::ios::out);

	if(nodes_file_stream.is_open() == true){
		for(int i = 0; i < quality_array.size(); i++)
			nodes_file_stream << quality_array[i] << std::endl;
		nodes_file_stream.close();
	}
	else
		std::cout << "File cannot be opened. " << std::endl;
}

bool compareQueueElement(queue_element e1, queue_element e2){
	return(e1.second < e2.second);
}

void ftkVesselTracer::GetBestTrace(std::vector<queue_element>& queue, queue_element& top){

	std::sort(queue.begin(), queue.end(), compareQueueElement);
	top.first = queue[0].first;
	top.second = queue[0].second;

	queue[0].second = this->allParams.nodeDetectionParams.infTraceQuality;
}

int ftkVesselTracer::TraceHitTest(Node node){

	int hit = 0;

	ImageType3D::IndexType node_index;
	node_index[0] = node.x; node_index[1] = node.y; node_index[2] = node.z;
	if(this->normalizedInputData->GetLargestPossibleRegion().IsInside(node_index) == false){
		hit = 1;
		return hit;
	}
	
	double norm = 0;
	bool parent_node_hit = false, common_parents = false;
	std::vector<double> parent_intersection(2 * node.parentIDLength, 0);
	
	for(int i = 0; i < this->allNodes.size(); i++){		
		
		parent_node_hit = false;
		for(int j = 0; j < node.parentID.size(); j++){
			if(i == node.parentID[j]){
				parent_node_hit = true;
				break;
			}
		}
		if(parent_node_hit == true)
			continue;

		norm = Node::ComputeNorm(Node(this->allNodes[i].x - node.x, this->allNodes[i].y - node.y, this->allNodes[i].z - node.z, 0));
		if(norm < this->allParams.nodeDetectionParams.distanceThresholdSecondary * (this->allNodes[i].scale + node.scale)){

			parent_intersection = std::vector<double>(2 * node.parentIDLength, 0);
			//std::set_intersection(this->allNodes[i].parentID.begin(), this->allNodes[i].parentID.end(), node.parentID.begin(), node.parentID.end(), parent_intersection.begin());
			//std::set_intersection(this->allNodes[i].parentID.begin(), this->allNodes[i].parentID.begin()+this->allNodes[i].parentID.size(), node.parentID.begin(), node.parentID.begin()+node.parentID.size(), parent_intersection.begin());
			//std::set_intersection(this->allNodes[i].parentID, this->allNodes[i].parentID+this->allNodes[i].parentID.size(), node.parentID, node.parentID+node.parentID.size(), parent_intersection.begin());
			
			common_parents = false;
			for(int j = 0; j < this->allNodes[i].parentID.size(); j++){
				for(int k = 0; k < node.parentID.size(); k++){
					if(this->allNodes[i].parentID[j] == node.parentID[k]){
						common_parents = true;
						break;
					}
				}
			}

			/*common_parents = false;
			for(int k = 0; k < parent_intersection.size(); k++){
				if(k != 0){
					common_parents = true;
					break;
				}
			}*/

			if(common_parents == true)
				continue;
			
			hit = 1;
			if(norm < std::min(this->allNodes[i].scale, node.scale))
				hit = 2;

			return hit;
		}
	}
	return hit;
}

double ftkVesselTracer::computeTraceQuality(Node& node){

	double cost = 0.0; 
	if(node.likelihood <= 0.0)
		cost = this->allParams.nodeDetectionParams.traceQualityThreshold; //.maxTraceCost;
	else
		cost = -1 * std::log(node.likelihood);

	int parent1 = node.parentID[0], parent2 = 0;
	std::vector<double> pos_diff(3, 0), pos_diff_last(3, 0);
	pos_diff[0] = this->allNodes[parent1].x - node.x;
	pos_diff[1] = this->allNodes[parent1].y - node.y;
	pos_diff[2] = this->allNodes[parent1].z - node.z;
	double pos_diff_norm = Node::ComputeNorm(Node(pos_diff[0], pos_diff[1], pos_diff[2], 0));
	pos_diff[0] = pos_diff[0] / pos_diff_norm;
	pos_diff[1] = pos_diff[1] / pos_diff_norm;
	pos_diff[2] = pos_diff[2] / pos_diff_norm;
	
	double current_curvature = 0;
	int count = 1;
	std::vector<double> arr_curvature;
	for(int i = 1; i < node.parentID.size(); i++){
		parent2 = node.parentID[i];

		if(parent2 > -1){	
			pos_diff_last = pos_diff;
			pos_diff[0] = this->allNodes[parent2].x - this->allNodes[parent1].x;
			pos_diff[1] = this->allNodes[parent2].y - this->allNodes[parent1].y;
			pos_diff[2] = this->allNodes[parent2].z - this->allNodes[parent1].z;
			pos_diff_norm = Node::ComputeNorm(Node(pos_diff[0], pos_diff[1], pos_diff[2], 0));
			pos_diff[0] = pos_diff[0] / pos_diff_norm;
			pos_diff[1] = pos_diff[1] / pos_diff_norm;
			pos_diff[2] = pos_diff[2] / pos_diff_norm;
			
			current_curvature = pos_diff[0] * pos_diff_last[0] + pos_diff[1] * pos_diff_last[1] + pos_diff[2] * pos_diff_last[2];
		}
		else
			current_curvature = 1.0;

		arr_curvature.push_back(current_curvature);
		
		// Cost function is different than one given in thesis?
		if(this->allNodes[parent1].likelihood <= 0.0)
			cost = this->allParams.nodeDetectionParams.traceQualityThreshold; //.maxTraceCost; //.traceQualityThreshold;
		else
			cost = cost - std::log(this->allNodes[parent1].likelihood);
		
		count = count + 1;

		if(parent2 == -1)
			break;
		parent1 = parent2;
	}

	double mean_curvature = std::accumulate(arr_curvature.begin(), arr_curvature.end(), 0.0);
	mean_curvature = mean_curvature / (double)arr_curvature.size();
	double std_curvature = 0.0;
	for(int i = 0; i < arr_curvature.size(); i++)
		std_curvature = std_curvature + std::pow(arr_curvature[i] - mean_curvature, 2);
	std_curvature = std::sqrt(std_curvature);

	node.traceQuality = std_curvature + (cost + this->allParams.nodeDetectionParams.traceLengthCost)/(double)count;

	return node.traceQuality;
}

void ftkVesselTracer::ComputeSecondaryNodeDirections(Node& node, std::vector<double>& dir_hist){
	
	// FUNCTION CALL
	node.nHoodSecondaryMultiplier = 2.0;

	node.nHoodScaleSecondary = node.nHoodSecondaryMultiplier * node.scale;
	double nHoodScale = node.nHoodScaleSecondary;
	int sMaxBin = this->allParams.oriBin.indexLength;
	nHoodScale = std::min(nHoodScale, (double)sMaxBin);

	ImageType3D::IndexType nhood_start_index, nhood_end_index;
		
	nhood_start_index[0] = floor(node.x - nHoodScale + 0.5); 
	nhood_start_index[1] = floor(node.y - nHoodScale + 0.5);
	nhood_start_index[2] = floor(node.z - nHoodScale + 0.5);

	nhood_end_index[0] = floor(node.x + nHoodScale + 0.5);
	nhood_end_index[1] = floor(node.y + nHoodScale + 0.5);
	nhood_end_index[2] = floor(node.z + nHoodScale + 0.5);
	
	//ImageType3D::IndexType current_index, current_index_1;
	int in_volume_count = 0, out_volume_count = 0; //, hist_bin = 0, max_bin_element = 0, min_bin_element = 0;
	//PixelType x_normalized = 0, y_normalized = 0, z_normalized = 0, current_pixel = 0;
	//double radial_dist = 0, hist_val_region_based = 0;
	//std::vector<int> bin_index(3, 0);
	std::vector<double> sph_hist_region_based(this->allParams.oriBin.angleCount, 0);
	std::vector<double> sph_hist_count(this->allParams.oriBin.angleCount, 0);

	//#pragma omp parallel for
	for(int j = nhood_start_index[1]; j <= nhood_end_index[1]; j++){
		for(int i = nhood_start_index[0]; i <= nhood_end_index[0]; i++){
			for(int k = nhood_start_index[2]; k <= nhood_end_index[2]; k++){
			//for(int j = nhood_start_index[1]; j <= nhood_end_index[1]; j++){
				
				ImageType3D::IndexType current_index, current_index_1;
				current_index[0] = i; current_index[1] = j; current_index[2] = k;
				//current_index_1[0] = j; current_index_1[1] = i; current_index_1[2] = k;
				
				if(this->normalizedInputData->GetLargestPossibleRegion().IsInside(current_index) == true){
					in_volume_count++;
				
					PixelType x_normalized = (i - node.x)/node.scale;
					PixelType y_normalized = (j - node.y)/node.scale;
					PixelType z_normalized = (k - node.z)/node.scale;
					double radial_dist = std::sqrt((x_normalized * x_normalized) + (y_normalized * y_normalized) + (z_normalized * z_normalized));

					//if(radial_dist > 1.0 && radial_dist < node.nHoodSecondaryMultiplier){
					if(radial_dist > 1.0 && radial_dist < node.nHoodSecondaryMultiplier	){
						
						std::vector<int> bin_index(3, 0);
						bin_index[0] = floor(i - node.x + 0.5); 
						bin_index[1] = floor(j - node.y + 0.5); 
						bin_index[2] = floor(k - node.z + 0.5);
						int max_bin_element = *std::max_element(bin_index.begin(), bin_index.end());
						int min_bin_element = *std::min_element(bin_index.begin(), bin_index.end());
						
						//#pragma omp critical
						if(max_bin_element < sMaxBin && min_bin_element > -1*sMaxBin){
							
							//int hist_bin1 = this->allParams.oriBin.binIndexVec[bin_index[0] + sMaxBin + 1][bin_index[1] + sMaxBin + 1][bin_index[2] + sMaxBin + 1];
							//int hist_bin4 = this->allParams.oriBin.binIndexVec[bin_index[0] + sMaxBin][bin_index[1] + sMaxBin][bin_index[2] + sMaxBin];

							//int hist_bin2 = this->allParams.oriBin.binIndexVec[bin_index[2] + sMaxBin][bin_index[1] + sMaxBin][bin_index[0] + sMaxBin];
							
							
							int hist_bin = this->allParams.oriBin.binIndexVec[bin_index[2] + sMaxBin][bin_index[0] + sMaxBin][bin_index[1] + sMaxBin];
							
							
							//int hist_bin3 = this->allParams.oriBin.binIndexVec[bin_index[0] + sMaxBin][bin_index[1] + sMaxBin][bin_index[2] + sMaxBin];
							//hist_bin = this->allParams.oriBin.binIndexVec[bin_index[1] + sMaxBin][bin_index[0] + sMaxBin][bin_index[2] + sMaxBin];
							
							PixelType current_pixel = this->normalizedInputData->GetPixel(current_index);
							//float current_pixel_1 = this->normalizedInputData->GetPixel(current_index_1);
							
							double hist_val_region_based = std::abs(current_pixel - node.meanBackgroundIntensity) - std::abs(current_pixel - node.meanForegroundIntensity);
							
							sph_hist_region_based[hist_bin] = sph_hist_region_based[hist_bin] + hist_val_region_based;
							sph_hist_count[hist_bin]++;
						}
					}
				}
				else
					out_volume_count++;
			}
		}
	}

	// Print the histogram
	//for(int i = 0; i < sph_hist_region_based.size(); i++)
	//	std::cout << i << ", " << sph_hist_region_based[i] << std::endl;

	ofstream nodes_file_stream;
	nodes_file_stream.open("SphHist_Cpp.txt", std::ios::out);

	if(nodes_file_stream.is_open() == true){
		for(int i = 0; i < sph_hist_region_based.size(); i++)
			nodes_file_stream << sph_hist_region_based[i] << std::endl;
		nodes_file_stream.close();
	}

	//#pragma omp parallel for
	// Normalize sphrical histogram
	for(int i = 0; i < sph_hist_region_based.size(); i++){
		sph_hist_region_based[i] = sph_hist_region_based[i]/(sph_hist_count[i] + 0.001);		
		sph_hist_region_based[i] = std::max(this->allParams.oriBin.minSphHistCount, sph_hist_region_based[i]);
	}

	double hist_sum = std::accumulate(sph_hist_region_based.begin(), sph_hist_region_based.end(), 0.0);
	if(hist_sum == 0){
		node.dirX.clear(); node.dirY.clear(); node.dirZ.clear();
		node.sphHistRegionBased = sph_hist_region_based;
		return;
	}
	//#pragma omp parallel for
	for(int i = 0; i < sph_hist_region_based.size(); i++)
		sph_hist_region_based[i] = sph_hist_region_based[i]/hist_sum;

	// Obtain the smoothing prior
	std::vector<double> smooth_hist_1(this->allParams.oriBin.angleCount, 0), smooth_hist_2(this->allParams.oriBin.angleCount, 0);
	for(int i = 0; i < this->allParams.oriBin.angleCount; i++){
		smooth_hist_1[i] = dir_hist[i];
		for(int j = 0; j < this->allParams.oriBin.nLastIndicesOfInterest; j++)
			smooth_hist_1[i] = smooth_hist_1[i] + this->allParams.oriBin.histSmoothingFactor * dir_hist[this->allParams.oriBin.nbr[i][j]];
	}
	for(int i = 0; i < this->allParams.oriBin.angleCount; i++){
		smooth_hist_2[i] = smooth_hist_1[i];
		for(int j = 0; j < this->allParams.oriBin.nLastIndicesOfInterest; j++)
			smooth_hist_2[i] = smooth_hist_2[i] + this->allParams.oriBin.histSmoothingFactor * smooth_hist_1[this->allParams.oriBin.nbr[i][j]];
	}
	
	//#pragma omp parallel for
	for(int i = 0; i < sph_hist_region_based.size(); i++)
		sph_hist_region_based[i] = sph_hist_region_based[i] * smooth_hist_2[i];
	
	hist_sum = std::accumulate(sph_hist_region_based.begin(), sph_hist_region_based.end(), 0.0);
	if(hist_sum == 0){
		node.dirX.clear(); node.dirY.clear(); node.dirZ.clear();
		node.sphHistRegionBased = sph_hist_region_based;
		return;
	}
	//#pragma omp parallel for
	for(int i = 0; i < sph_hist_region_based.size(); i++)
		sph_hist_region_based[i] = sph_hist_region_based[i]/hist_sum;

	// Print the histogram
	for(int i = 0; i < sph_hist_region_based.size(); i++)
		std::cout << i << ", " << sph_hist_region_based[i] << std::endl;
	
	ofstream nodes_file_stream1;
	nodes_file_stream1.open("SphHist_Cpp_smoothed.txt", std::ios::out);

	if(nodes_file_stream1.is_open() == true){
		for(int i = 0; i < sph_hist_region_based.size(); i++)
			nodes_file_stream1 << sph_hist_region_based[i] << std::endl;
		nodes_file_stream1.close();
	}
	
	// Finding modes of the histogram
	bool mode_flag = false;
	std::vector<std::vector<double> > mode_bins;
	std::vector<double> mode_hist_val;
	std::vector<double> a_bin(3, 0);
	int index1 = 0, index2 = 0;
	// Find out the peaks of the histogram and the corresponding bin locations
	for(int i = 0; i < this->allParams.oriBin.angleCount; i++){
		mode_flag = true;
		if(sph_hist_region_based[i] == 0)
			continue;
		
		for(int j = 0; j < this->allParams.oriBin.nLastIndicesOfInterest; j++){
			
			//CHECK THE ACCESSING OF ALL 2D ARRAYS (COMPARED TO MATLAB)
			index1 = this->allParams.oriBin.nbr[i][j];
			index2 = i;
			
			//if(sph_hist_region_based[this->allParams.oriBin.nbr[i][j]] > sph_hist_region_based[i]){
			if(sph_hist_region_based[index1] > sph_hist_region_based[index2]){
				double v1 = sph_hist_region_based[index1];
				double v2 = sph_hist_region_based[index2];
				mode_flag = false;
				break;
			}
		}

		if(mode_flag == true){	
			a_bin[0] = this->allParams.oriBin.binCenters[0][i];
			a_bin[1] = this->allParams.oriBin.binCenters[1][i];
			a_bin[2] = this->allParams.oriBin.binCenters[2][i];
			mode_bins.push_back(a_bin);
			mode_hist_val.push_back(sph_hist_region_based[i]);
		}
	}
	
	// Refine node directions: get the top 2 peaks separated by 60 deg, see if the 3rd also higher than
	// 0.25*(p1 + p2)
	std::vector<std::vector<double> > mode_bins_sorted;
	std::vector<double> mode_hist_val_sorted;
	std::multimap<double, std::vector<double> > bin_and_val_map;
	std::multimap<double, std::vector<double> >::iterator iter_map;
	for(int i = 0; i < mode_bins.size(); i++)
		bin_and_val_map.insert(std::pair<double, std::vector<double> >(mode_hist_val[i], mode_bins[i]));
	
	for(iter_map = bin_and_val_map.begin(); iter_map != bin_and_val_map.end(); iter_map++){
		mode_hist_val_sorted.push_back((*iter_map).first);
		mode_bins_sorted.push_back((*iter_map).second);
	}
	std::reverse(mode_hist_val_sorted.begin(), mode_hist_val_sorted.end());
	std::reverse(mode_bins_sorted.begin(), mode_bins_sorted.end());
	
	double x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0;
	std::vector<double> bin1(3, 0), bin2(3, 0);
	for(int i = 1; i < mode_bins_sorted.size(); i++){
		bin1 = mode_bins_sorted[i];
		for(int j = 0; j <= i - 1; j++){
			bin2 = mode_bins_sorted[j];
			if((bin1[0]*bin2[0] + bin1[1]*bin2[1] + bin1[2]*bin2[2]) > std::cos(this->allParams.nodeDetectionParams.maxBranchAngle*vnl_math::pi/180.0))
				mode_hist_val_sorted[i] = 0;
		}
	}

	std::vector<std::vector<double> > mode_bins_sorted_1;
	std::vector<double> mode_hist_val_sorted_1;
	std::multimap<double, std::vector<double> > bin_and_val_map_1;
	std::multimap<double, std::vector<double> >::iterator iter_map_1;
	for(int i = 0; i < mode_bins_sorted.size(); i++)
		bin_and_val_map_1.insert(std::pair<double, std::vector<double> >(mode_hist_val_sorted[i], mode_bins_sorted[i]));
	
	bin_and_val_map_1.erase(0);

	for(iter_map_1 = bin_and_val_map_1.begin(); iter_map_1 != bin_and_val_map_1.end(); iter_map_1++){
		mode_hist_val_sorted_1.push_back((*iter_map_1).first);
		mode_bins_sorted_1.push_back((*iter_map_1).second);
	}
	std::reverse(mode_hist_val_sorted_1.begin(), mode_hist_val_sorted_1.end());
	std::reverse(mode_bins_sorted_1.begin(), mode_bins_sorted_1.end());
	
	double branching_th = 0;
	std::vector<double> dirX, dirY, dirZ;
	if(mode_bins_sorted_1.size() > 2){
		branching_th = this->allParams.nodeDetectionParams.branchingThreshold * (mode_hist_val_sorted_1[0] + mode_hist_val_sorted_1[1]);

		if(mode_hist_val_sorted_1[2] > branching_th){
			for(int i = 0; i < 3; i++){
				dirX.push_back(mode_bins_sorted_1[i][0]);
				dirY.push_back(mode_bins_sorted_1[i][1]);
				dirZ.push_back(mode_bins_sorted_1[i][2]);
			}
		}
		else{
			for(int i = 0; i < 2; i++){
				dirX.push_back(mode_bins_sorted_1[i][0]);
				dirY.push_back(mode_bins_sorted_1[i][1]);
				dirZ.push_back(mode_bins_sorted_1[i][2]);
			}
		}
	}
	else{
		for(int i = 0; i < mode_bins_sorted_1.size(); i++){
			dirX.push_back(mode_bins_sorted_1[i][0]);
			dirY.push_back(mode_bins_sorted_1[i][1]);
			dirZ.push_back(mode_bins_sorted_1[i][2]);
		}
	}
	if(dirX.size() < 2){
		dirX.push_back(-1 * dirX[0]);
		dirY.push_back(-1 * dirY[0]);
		dirZ.push_back(-1 * dirZ[0]);
	}

	node.dirX = dirX;
	node.dirY = dirY;
	node.dirZ = dirZ;
	node.sphHistRegionBased = sph_hist_region_based;

	// Print the detected maximal nodes
	//for(int i = 0; i < node.dirX.size(); i++){
	//	std::cout << node.dirX[i] << ", " << node.dirY[i] << ", " << node.dirZ[i] << std::endl;
	//}
}

void ftkVesselTracer::writeNodesToFile(std::vector<Node> nodes, std::string& file_path){

	ofstream nodes_file_stream;
	nodes_file_stream.open(file_path.c_str(), std::ios::out);

	if(nodes_file_stream.is_open() == true){

		for(int i = 0; i < nodes.size(); i++){

			nodes_file_stream << nodes[i].x << "," << nodes[i].y << "," << nodes[i].z << "," << nodes[i].scale << "," << nodes[i].likelihood << "," 
				<< nodes[i].meanForegroundIntensity << "," << nodes[i].meanBackgroundIntensity << ",";
			for(int j = 0; j < nodes[i].parentID.size(); j++)
				nodes_file_stream << nodes[i].parentID[j] << ",";
			nodes_file_stream << std::endl;
		}
		nodes_file_stream.close();
	}
	else
		std::cout << "Unable to open file for writing nodes: " << file_path.c_str() << std::endl;
}

void ftkVesselTracer::ReadNodesFromTextFile(const std::string& file_name){

	if(this->allNodes.size() != 0)
		this->allNodes.clear();

	ifstream nodes_file_stream;
	nodes_file_stream.open(file_name.c_str(), std::ios::in);
	
	std::string str, str1;
	std::vector<std::string> node_str;
	if(nodes_file_stream.is_open() == true){

		while(nodes_file_stream.good() == true){
			
			std::string str;
			if(getline(nodes_file_stream, str) == false)
				break;
			//std::cout << str.c_str() << std::endl;

			std::istringstream str_stream(str); 
			while(str_stream.good() == true){
				
				if(getline(str_stream, str1, ',') == false)
					break;
				node_str.push_back(str1);
			}

			Node n1;
			n1.x = atof(node_str[0].c_str());
			n1.y = atof(node_str[1].c_str());
			n1.z = atof(node_str[2].c_str());
			n1.scale = atof(node_str[3].c_str());
			n1.likelihood = atof(node_str[4].c_str());
			n1.meanForegroundIntensity = atof(node_str[5].c_str());
			n1.meanBackgroundIntensity = atof(node_str[6].c_str());
			
			// Works only for parentID size = 4
			n1.parentIDLength = 4;
			n1.parentID = std::vector<double>(n1.parentIDLength, 0);
			n1.parentID[0] = atoi(node_str[7].c_str());
			n1.parentID[1] = atoi(node_str[8].c_str());
			n1.parentID[2] = atoi(node_str[9].c_str());
			n1.parentID[3] = atoi(node_str[10].c_str());
			
			this->allNodes.push_back(n1);	

			node_str.clear();
		}
		//if(!this->allNodes.empty())
			//this->VisualizeNodesWithData3D(this->allNodes, false);
	}
	else
		std::cout << "Unable to open the nodes file: " << file_name.c_str() << std::endl;
}

void Node::ComputeDistanceBetweenNodes(Node n1, Node n2, Node &nd){

	nd.x = n1.x - n2.x;
	nd.y = n1.y - n2.y;
	nd.z = n1.z - n2.z;
}

void Node::NormalizeNode(double norm_val){
	
	this->x = this->x / norm_val;
	this->y = this->y / norm_val;
	this->z = this->z / norm_val;
}

void Node::InvertNodeDir(void){

	this->x = -1 * this->x;
	this->y = -1 * this->y;
	this->z = -1 * this->z;
}

double Node::DotProduct(Node n1, Node n2){

	return(n1.x*n2.x + n1.y*n2.y + n1.z*n2.z);
}

AffinityEdge::AffinityEdge(void){

	this->from = 0;
	this->to = 0;
	this->weight = 0.0;
}

AffinityEdge::AffinityEdge(int from_node_index, int to_node_index, double edge_weight){

	this->from = from_node_index;
	this->to = to_node_index;
	this->weight = edge_weight;
}

bool CompareEdges(AffinityEdge e1, AffinityEdge e2){
	return(e1.weight > e2.weight);
}

Tree::Tree(void){

	this->ID = 0;
	this->start = 0;
	this->NNodes = 0;
}

Tree::Tree(int ID, int start, int NNodes){

	this->ID = ID;
	this->start = start;
	this->NNodes = NNodes;
}

void ftkVesselTracer::CreateMinimumSpanningForest(){
	
	this->CreateAffinityGraph();
	//this->ComputeMinimumSpanningTreeBoost();
	this->ComputeMinimumSpanningForestWithLoopDetection();
}

void ftkVesselTracer::CreateAffinityGraph(void){

	//Node dir;
	//double norm_dist = 0.0, existingNodeDist = 0.0;
	//int angleBinLinear = 0;
	
	//allocate memory for the binned edges in every node
	for(int i = 0; i < this->allNodes.size(); i++)
		this->allNodes[i].connectedNodesBinned.resize(2 * this->allParams.graphAndMSTParams.NBinsAffinity, this->allParams.graphAndMSTParams.maxEdgeWeight+1);
	
	for(int i = 0; i < this->allNodes.size(); i++){
		
		#pragma omp parallel for
		for(int j = i + 1; j < this->allNodes.size(); j++){
			
			Node dir;
			Node::ComputeDistanceBetweenNodes(this->allNodes[j], this->allNodes[i], dir);
			double norm_dist = Node::ComputeNorm(dir);
			dir.NormalizeNode(norm_dist);

			if(norm_dist < (1.0) * this->allParams.graphAndMSTParams.affinityRadThresh * std::max(this->allNodes[j].scale, this->allNodes[i].scale)){
				
				int angleBinLinear = this->ComputeAffinityBin(dir);
				double existingNodeDist = this->allNodes[i].connectedNodesBinned[2*(angleBinLinear - 1) + 1];

				if(existingNodeDist > norm_dist){
					this->allNodes[i].connectedNodesBinned[2*(angleBinLinear - 1)] = j;
					this->allNodes[i].connectedNodesBinned[2*(angleBinLinear - 1) + 1] = norm_dist;
				}
				
				dir.InvertNodeDir();

				angleBinLinear = this->ComputeAffinityBin(dir);
				existingNodeDist = this->allNodes[j].connectedNodesBinned[2*(angleBinLinear - 1) + 1];
				
				if(existingNodeDist > norm_dist){
					this->allNodes[j].connectedNodesBinned[2*(angleBinLinear - 1)] = i;
					this->allNodes[j].connectedNodesBinned[2*(angleBinLinear - 1) + 1] = norm_dist;
				}				
			}
		}
	}
	
	Node dir1, dir;
	double norm_dist1 = 0.0, edge_weight = 0.0;
	double norm_dist = 0.0, existingNodeDist = 0.0;
	for(int i = 0; i < this->allNodes.size(); i++){
		for(int j = 0; j < 2*this->allParams.graphAndMSTParams.NBinsAffinity; j = j+2){
			if(this->allNodes[i].connectedNodesBinned[j+1] > this->allParams.graphAndMSTParams.maxEdgeWeight)
				continue;

			Node::ComputeDistanceBetweenNodes(this->allNodes[i], this->allNodes[this->allNodes[i].connectedNodesBinned[j]], dir);
			norm_dist = Node::ComputeNorm(dir);
			dir.NormalizeNode(norm_dist);

			for(int k = j+2; k < 2*this->allParams.graphAndMSTParams.NBinsAffinity; k = k+2){
				if(this->allNodes[i].connectedNodesBinned[k+1] > this->allParams.graphAndMSTParams.maxEdgeWeight)
					continue;
				
				Node::ComputeDistanceBetweenNodes(this->allNodes[i], this->allNodes[this->allNodes[i].connectedNodesBinned[k]], dir1);
				norm_dist1 = Node::ComputeNorm(dir1);
				dir1.NormalizeNode(norm_dist1);
				
				if(Node::DotProduct(dir, dir1) > cos(this->allParams.graphAndMSTParams.minBranchAngle)){
					if(this->allNodes[i].connectedNodesBinned[j+1] > this->allNodes[i].connectedNodesBinned[k+1])
						this->allNodes[i].connectedNodesBinned[j+1] = this->allParams.graphAndMSTParams.maxEdgeWeight + 1;
					else
						this->allNodes[i].connectedNodesBinned[k+1] = this->allParams.graphAndMSTParams.maxEdgeWeight + 1;
				}
			}
		}

		for(int j = 0; j < 2*this->allParams.graphAndMSTParams.NBinsAffinity; j = j+2){
			if(this->allNodes[i].connectedNodesBinned[j+1] < this->allParams.graphAndMSTParams.maxEdgeWeight + 1){
				
				// Edge weights are not similar to those mentioned in thesis
				// This edge weight is inverse of Amit's Matlab code - Higher weight would mean more likelihood of an edge being present
				//edge_weight = this->allNodes[i].connectedNodesBinned[j+1] / std::min(this->allNodes[i].likelihood, this->allNodes[this->allNodes[i].connectedNodesBinned[j]].likelihood);
				edge_weight = std::min(this->allNodes[i].likelihood, this->allNodes[this->allNodes[i].connectedNodesBinned[j]].likelihood)/this->allNodes[i].connectedNodesBinned[j+1];
				
				this->allNodes[i].connectedNodesAffinity.push_back(std::pair<int, double>(this->allNodes[i].connectedNodesBinned[j], edge_weight));

				this->edges.push_back(AffinityEdge(i, this->allNodes[i].connectedNodesBinned[j], edge_weight));
			}
		}
	}

	//this->VisualizeAffinityGraph();



	//Toy graph
	/*this->edges.clear();
	this->edges.push_back(AffinityEdge(1, 7, 5));
	this->edges.push_back(AffinityEdge(7, 1, 5));
	this->edges.push_back(AffinityEdge(7, 6, 2));
	this->edges.push_back(AffinityEdge(6, 7, 2));
	this->edges.push_back(AffinityEdge(7, 4, 10));
	this->edges.push_back(AffinityEdge(4, 7, 10));
	this->edges.push_back(AffinityEdge(4, 6, 7));
	this->edges.push_back(AffinityEdge(6, 4, 7));
	this->edges.push_back(AffinityEdge(4, 5, 5));
	this->edges.push_back(AffinityEdge(5, 4, 5));
	this->edges.push_back(AffinityEdge(10, 2, 4));
	this->edges.push_back(AffinityEdge(2, 10, 4));
	this->edges.push_back(AffinityEdge(2, 3, 5));
	this->edges.push_back(AffinityEdge(3, 2, 5));
	this->edges.push_back(AffinityEdge(2, 9, 12));
	this->edges.push_back(AffinityEdge(9, 2, 12));
	this->edges.push_back(AffinityEdge(2, 8, 3));
	this->edges.push_back(AffinityEdge(8, 2, 3));
	this->edges.push_back(AffinityEdge(3, 9, 5));
	this->edges.push_back(AffinityEdge(9, 3, 5));
	this->edges.push_back(AffinityEdge(9, 8, 2));
	this->edges.push_back(AffinityEdge(8, 9, 2));
	this->edges.push_back(AffinityEdge(10, 4, 30));
	this->edges.push_back(AffinityEdge(4, 10, 30));*/


	

	// Sort the edges in descending order of weights (since higher weight is better)
	std::sort(this->edges.begin(), this->edges.end(), CompareEdges);

	// Initialize all forest nodes
	this->allForestNodes = std::vector<Node>(this->allNodes.size());
	for(int i = 0; i < this->allNodes.size(); i++){
		this->allForestNodes[i].x = this->allNodes[i].x;
		this->allForestNodes[i].y = this->allNodes[i].y;
		this->allForestNodes[i].z = this->allNodes[i].z;
		this->allForestNodes[i].likelihood = this->allNodes[i].likelihood;

		this->allForestNodes[i].ID = this->allNodes[i].ID;
		this->allForestNodes[i].branchIDs = std::vector<int>(this->allParams.graphAndMSTParams.maxNBranches, -1);
		this->allForestNodes[i].NBranches = this->allNodes[i].NBranches;
	}

	// Used no longer
	//this->allNodes.clear();
}

int ftkVesselTracer::ComputeAffinityBin(Node dir){

	double phi = (atan2(dir.y, dir.x) / vnl_math::pi) + 1;
	double theta = acos(dir.z) / vnl_math::pi; 

	phi = floor(4*phi + 0.5);
	theta = floor(4*theta + 0.5);

	if(phi == 0)
		phi = 8;
	if(theta == 0)
		theta = 4;

	int bin = phi + 8*(theta - 1);
	
	return bin;
}

void ftkVesselTracer::VisualizeAffinityGraph(void){

	ITKToVTKConnectorType::Pointer ITK_to_VTK_connector = ITKToVTKConnectorType::New();

	ITK_to_VTK_connector->SetInput(this->inputDataForRendering);
	ITK_to_VTK_connector->Update();

	vtkSmartPointer<vtkImageData> vtk_image = ITK_to_VTK_connector->GetOutput();

	// Testing vtk image
	//vtk_image->PrintSelf(std::cout, vtkIndent(0));

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(1.0, 1.0, 1.0);
	
	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	renderer->SetActiveCamera(camera);

	vtkSmartPointer<vtkRenderWindow> render_window = vtkSmartPointer<vtkRenderWindow>::New();
	render_window->AddRenderer(renderer);
	
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	render_window_interactor->SetRenderWindow(render_window);
	
	vtkSmartPointer<vtkPiecewiseFunction> opacity_transfer_function = vtkSmartPointer<vtkPiecewiseFunction>::New();
	opacity_transfer_function->AddPoint(2, 0.0);
	opacity_transfer_function->AddPoint(10, 0.1);
	
	vtkSmartPointer<vtkColorTransferFunction> color_transfer_function = vtkSmartPointer<vtkColorTransferFunction>::New();
	color_transfer_function->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	color_transfer_function->AddRGBPoint((10 * 255)/100, 1, 1, 1); //blue
	color_transfer_function->AddRGBPoint((45 * 255)/100, 0, .01, 0); //green
	color_transfer_function->AddRGBPoint((150 * 255)/100, .01, 0, 0); //red
	
	vtkSmartPointer<vtkVolumeProperty> volume_property = vtkSmartPointer<vtkVolumeProperty>::New();
	volume_property->SetColor(color_transfer_function);
	volume_property->SetScalarOpacity(opacity_transfer_function);
	volume_property->ShadeOn();
	volume_property->SetInterpolationTypeToLinear();
	
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volume_mapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	
	volume_mapper->SetInput(vtk_image);
	volume_mapper->SetBlendModeToComposite();
	
	vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(volume_mapper); 
	volume->SetProperty(volume_property);
	volume->SetPosition(0, 0, 0);
	volume->SetPickable(0);
	volume->Update();


	// Preparing vtkPolyData to show the graph
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	
	int point_counter = -1;
	for(int i = 0; i < this->allNodes.size(); i++){
		points->InsertNextPoint(this->allNodes[i].x, this->allNodes[i].y, this->allNodes[i].z);
		point_counter++;

		for(int j = 0; j < this->allNodes[i].connectedNodesAffinity.size(); j++){
			points->InsertNextPoint(this->allNodes[this->allNodes[i].connectedNodesAffinity[j].first].x, 
				this->allNodes[this->allNodes[i].connectedNodesAffinity[j].first].y, this->allNodes[this->allNodes[i].connectedNodesAffinity[j].first].z);			
			point_counter++;
		}
		for(int j = 0; j < this->allNodes[i].connectedNodesAffinity.size(); j++){
			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
			line->GetPointIds()->SetId(0, point_counter - this->allNodes[i].connectedNodesAffinity.size());
			line->GetPointIds()->SetId(1, point_counter - j);
			lines->InsertNextCell(line);
		}
	}
	
	unsigned char green[3] = {0, 255, 0};
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->InsertNextTupleValue(green);
	
	vtkSmartPointer<vtkPolyData> poly_data = vtkSmartPointer<vtkPolyData>::New();
	poly_data->SetPoints(points);
	poly_data->SetLines(lines);
	//poly_data->GetCellData()->SetScalars(colors);
	
	vtkSmartPointer<vtkPolyDataMapper> poly_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	poly_mapper->SetInput(poly_data);

	vtkSmartPointer<vtkActor> poly_actor = vtkSmartPointer<vtkActor>::New();
	poly_actor->SetMapper(poly_mapper);
	poly_actor->GetProperty()->SetLineWidth(3);
	poly_actor->GetProperty()->SetColor(0, 255, 0);
	
	renderer->AddActor(poly_actor);
	renderer->AddVolume(volume);

	renderer->ResetCamera();
	
	render_window_interactor->Initialize();
	render_window->Render();
	render_window_interactor->Start();
}

void ftkVesselTracer::VisualizeMinimumSpanningForest(void){
	
	ITKToVTKConnectorType::Pointer ITK_to_VTK_connector = ITKToVTKConnectorType::New();

	ITK_to_VTK_connector->SetInput(this->inputDataForRendering);
	ITK_to_VTK_connector->Update();

	vtkSmartPointer<vtkImageData> vtk_image = ITK_to_VTK_connector->GetOutput();

	// Testing vtk image
	//vtk_image->PrintSelf(std::cout, vtkIndent(0));

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(1.0, 1.0, 1.0);
	
	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	renderer->SetActiveCamera(camera);

	vtkSmartPointer<vtkRenderWindow> render_window = vtkSmartPointer<vtkRenderWindow>::New();
	render_window->AddRenderer(renderer);
	
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	render_window_interactor->SetRenderWindow(render_window);
	
	vtkSmartPointer<vtkPiecewiseFunction> opacity_transfer_function = vtkSmartPointer<vtkPiecewiseFunction>::New();
	opacity_transfer_function->AddPoint(2, 0.0);
	opacity_transfer_function->AddPoint(10, 0.1);
	
	vtkSmartPointer<vtkColorTransferFunction> color_transfer_function = vtkSmartPointer<vtkColorTransferFunction>::New();
	color_transfer_function->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	color_transfer_function->AddRGBPoint((10 * 255)/100, 1, 1, 1); //blue
	color_transfer_function->AddRGBPoint((45 * 255)/100, 0, .01, 0); //green
	color_transfer_function->AddRGBPoint((150 * 255)/100, .01, 0, 0); //red
	
	vtkSmartPointer<vtkVolumeProperty> volume_property = vtkSmartPointer<vtkVolumeProperty>::New();
	volume_property->SetColor(color_transfer_function);
	volume_property->SetScalarOpacity(opacity_transfer_function);
	volume_property->ShadeOn();
	volume_property->SetInterpolationTypeToLinear();
	
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volume_mapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	
	volume_mapper->SetInput(vtk_image);
	volume_mapper->SetBlendModeToComposite();
	
	vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(volume_mapper); 
	volume->SetProperty(volume_property);
	volume->SetPosition(0, 0, 0);
	volume->SetPickable(0);
	volume->Update();


	// Preparing vtkPolyData to show the graph
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	
	int point_counter = -1;
	for(int i = 0; i < this->allForestNodes.size(); i++){
		points->InsertNextPoint(this->allForestNodes[i].x, this->allForestNodes[i].y, this->allForestNodes[i].z);
		point_counter++;

		for(int j = 0; j < this->allForestNodes[i].branchIDs.size(); j++){
			points->InsertNextPoint(this->allForestNodes[this->allForestNodes[i].branchIDs[j]].x, 
				this->allForestNodes[this->allForestNodes[i].branchIDs[j]].y, this->allForestNodes[this->allForestNodes[i].branchIDs[j]].z);			
			point_counter++;
		}
		for(int j = 0; j < this->allForestNodes[i].branchIDs.size(); j++){
			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
			line->GetPointIds()->SetId(0, point_counter - this->allForestNodes[i].branchIDs.size());
			line->GetPointIds()->SetId(1, point_counter - j);
			lines->InsertNextCell(line);
		}
	}
	
	unsigned char green[3] = {0, 255, 0};
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->InsertNextTupleValue(green);
	
	vtkSmartPointer<vtkPolyData> poly_data = vtkSmartPointer<vtkPolyData>::New();
	poly_data->SetPoints(points);
	poly_data->SetLines(lines);
	//poly_data->GetCellData()->SetScalars(colors);
	
	vtkSmartPointer<vtkPolyDataMapper> poly_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	poly_mapper->SetInput(poly_data);

	vtkSmartPointer<vtkActor> poly_actor = vtkSmartPointer<vtkActor>::New();
	poly_actor->SetMapper(poly_mapper);
	poly_actor->GetProperty()->SetLineWidth(3);
	poly_actor->GetProperty()->SetColor(0, 255, 0);
	
	renderer->AddActor(poly_actor);
	renderer->AddVolume(volume);

	renderer->ResetCamera();
	
	render_window_interactor->Initialize();
	render_window->Render();
	render_window_interactor->Start() ;

}
void ftkVesselTracer::ComputeMinimumSpanningTreeBoost(void){

	//vtkSmartPointer<vtkMutableUndirectedGraph> vtk_graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
}

void ftkVesselTracer::ComputeMinimumSpanningForestWithLoopDetection(void){

	
	int from = 0, to = 0, label = -1, oldID = 0, newID = 0, loop_count = 0;
	for(int i = 0; i < this->edges.size(); i++){		
		from = this->edges[i].from;
		to = this->edges[i].to;

		if(this->allForestNodes[from].ID  == -1 && this->allForestNodes[to].ID == -1){
			
			label++;
			
			this->forest.push_back(Tree(label, from, 2));

			this->allForestNodes[from].ID = label;
			this->allForestNodes[to].ID = label;

			this->allForestNodes[from].NBranches++;
			this->allForestNodes[from].branchIDs[this->allForestNodes[from].NBranches - 1] = to;

			this->allForestNodes[to].NBranches++;
			this->allForestNodes[to].branchIDs[this->allForestNodes[to].NBranches - 1] = from;
		}
		else if(this->allForestNodes[from].ID == -1 && this->allForestNodes[to].ID != -1){
			
			this->allForestNodes[from].ID = this->allForestNodes[to].ID;
			
			this->allForestNodes[from].NBranches++;
			this->allForestNodes[from].branchIDs[this->allForestNodes[from].NBranches - 1] = to;

			this->allForestNodes[to].NBranches++;
			this->allForestNodes[to].branchIDs[this->allForestNodes[to].NBranches - 1] = from;

			this->forest[this->allForestNodes[to].ID].NNodes++;
		}
		else if(this->allForestNodes[from].ID != -1 && this->allForestNodes[to].ID == -1){
			
			this->allForestNodes[to].ID = this->allForestNodes[from].ID;
			
			this->allForestNodes[from].NBranches++;
			this->allForestNodes[from].branchIDs[this->allForestNodes[from].NBranches - 1] = to;

			this->allForestNodes[to].NBranches++;
			this->allForestNodes[to].branchIDs[this->allForestNodes[to].NBranches - 1] = from;

			this->forest[this->allForestNodes[from].ID].NNodes++;
		}
		else{
			if(this->allForestNodes[to].ID < this->allForestNodes[from].ID){
				
				oldID = this->allForestNodes[from].ID;
				newID = this->allForestNodes[to].ID;
				
				this->RelabelForestNodes(oldID, newID);

				this->forest[newID].NNodes = this->forest[newID].NNodes + this->forest[oldID].NNodes;
				this->forest[oldID].NNodes = 0;

				this->allForestNodes[from].NBranches++;
				this->allForestNodes[from].branchIDs[this->allForestNodes[from].NBranches - 1] = to;

				this->allForestNodes[to].NBranches++;
				this->allForestNodes[to].branchIDs[this->allForestNodes[to].NBranches - 1] = from;
			}
			else if(this->allForestNodes[to].ID > this->allForestNodes[from].ID){

				oldID = this->allForestNodes[to].ID;
				newID = this->allForestNodes[from].ID;

				this->RelabelForestNodes(oldID, newID);

				this->forest[newID].NNodes = this->forest[newID].NNodes + this->forest[oldID].NNodes;
				this->forest[oldID].NNodes = 0;
				
				this->allForestNodes[from].NBranches++;
				this->allForestNodes[from].branchIDs[this->allForestNodes[from].NBranches - 1] = to;

				this->allForestNodes[to].NBranches++;
				this->allForestNodes[to].branchIDs[this->allForestNodes[to].NBranches - 1] = from;
			}
			else if(this->CheckNeighbors(from, to) == false){
				loop_count++;
				this->loops.push_back(AffinityEdge(from, to, 0.0));	
			}
		}
	}
	
	std::cout << "Forest created with " << this->edges.size() - loop_count << " edges. " << std::endl;

	//this->PrintForest();

	//Removing redundant branches
	for(int i = 0; i < this->allForestNodes.size(); i++){
		for(int j = 0; j < this->allForestNodes[i].branchIDs.size(); j++){
			if(this->allForestNodes[i].branchIDs[j] == -1){
				this->allForestNodes[i].branchIDs.erase(this->allForestNodes[i].branchIDs.begin() + j, this->allForestNodes[i].branchIDs.end());
				break;
			}
		}
	}

	this->VisualizeMinimumSpanningForest();

	/*int count = -1;
	for(int i = 0; i < label; i++){
		if(this->forest[i].NNodes > 0){
			count++;
			this->forest[count].NNodes = this->forest[i].NNodes;
			this->forest[count].start = this->forest[i].start;

			this->forest[count].ID = i;
		}
		else
			this->forest.erase(this->forest.begin() + i);
	}

	this->VisualizeMinimumSpanningForest();*/
}

void ftkVesselTracer::RelabelForestNodes(int oldID, int newID){
	
	std::vector<int> connected_nodes;
	this->GetTree(oldID, connected_nodes);

	if(connected_nodes.empty() == false){
		for(int i = 0; i < connected_nodes.size(); i++)
			this->allForestNodes[connected_nodes[i]].ID = newID;
	}
}

void ftkVesselTracer::GetTree(int ID2, std::vector<int>& connected_nodes){

	int root = this->forest[ID2].start, ID = this->forest[ID2].ID, currentNNodes = 0, totalNNodes = 0, childNodeID = 0;
	
	if(this->allForestNodes[root].ID != ID)
		return;
	
	connected_nodes.resize(this->allParams.graphAndMSTParams.maxTreeNodes + 100, 0);
	connected_nodes[currentNNodes] = root;

	while(currentNNodes <= totalNNodes){
		
		for(int i = 0; i < this->allForestNodes[connected_nodes[currentNNodes]].NBranches; i++){
			childNodeID = this->allForestNodes[connected_nodes[currentNNodes]].branchIDs[i];
			
			//if(this->allForestNodes[childNodeID].ID == ID && std::find(connected_nodes.begin(), connected_nodes.begin() + totalNNodes, childNodeID) == connected_nodes.end()){
			if(this->allForestNodes[childNodeID].ID == ID && ftkVesselTracer::IsInList(connected_nodes, totalNNodes, childNodeID) == false){
				totalNNodes++;
				connected_nodes[totalNNodes] = childNodeID;
			}
		}
		currentNNodes++;
	}

	connected_nodes.erase(connected_nodes.begin() + totalNNodes + 1, connected_nodes.end());

	if(this->forest[ID2].NNodes != connected_nodes.size()){
		std::cout << "Error in scanning tree ID = " << ID << " NNodes = " << this->forest[ID2].NNodes << " Connected nodes = ";
		for(int i = 0; i < connected_nodes.size(); i++)
			std::cout << " " << connected_nodes[i];
		std::cout << std::endl;
	}
}

bool ftkVesselTracer::CheckNeighbors(int ID1, int ID2){

	bool is_neighbor = false;
	
	for(int i = 0; i < this->allForestNodes[ID1].NBranches; i++){
		if(this->allForestNodes[ID1].branchIDs[i] == ID2)
			is_neighbor = true;
	}

	for(int i = 0; i < this->allForestNodes[ID2].NBranches; i++){
		if(this->allForestNodes[ID2].branchIDs[i] == ID1)
			is_neighbor = true;
	}
	return is_neighbor;
}

bool ftkVesselTracer::IsInList(std::vector<int> list, int last_item, int item_to_find){
	
	bool present = false;
	for(int i = 0; i <= last_item; i++){
		if(list[i] == item_to_find){
			present = true;
			break;
		}	
	}
	return present;
}

void ftkVesselTracer::PrintForest(void){
	
	std::vector<int> connected_nodes; 
	for(int i = 0; i < this->forest.size(); i++){
		std::cout << " Tree: " << i << " with " << this->forest[i].NNodes << " elements [";
		
		if(connected_nodes.empty() == false)
			connected_nodes.clear();
		
		this->GetTree(i, connected_nodes);

		if(connected_nodes.empty() == false){
			for(int j = 0; j < connected_nodes.size(); j++)
				std::cout << connected_nodes[j] << " ";
			std::cout << "] " << std::endl;
		}
		else
			std::cout << "NULL]" << std::endl;
	}
}