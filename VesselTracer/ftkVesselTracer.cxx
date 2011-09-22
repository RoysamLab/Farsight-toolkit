
#include "ftkVesselTracer.h"

ftkVesselTracer::ftkVesselTracer(){
}

ftkVesselTracer::ftkVesselTracer(std::string input_data_path, bool preprocess = false){

	if(preprocess){
		ImageType3D::Pointer tiff_data_ptr; //only required for preprocessing
		this->PreprocessData(input_data_path, tiff_data_ptr); // ask user	
	}

	// Load preprocessed data
	this->LoadPreprocessedData(input_data_path);

	// 3D renering of original and preprocessed data
	Common::RescaleDataForRendering(this->originalData, this->originalDataForRendering);
	Common::RescaleDataForRendering(this->inputData, this->inputDataForRendering);
	/*Common::RenderImage3D(this->originalDataForRendering);
	Common::RenderImage3D(this->inputDataForRendering);
	*/

	// MIP computation
	this->ComputeMaximumIntensityProjectionImage();
	this->RenderMIPImage();	
}
ftkVesselTracer::~ftkVesselTracer(){
}

int ftkVesselTracer::PreprocessData(std::string file_path, ImageType3D::Pointer& data_ptr){

	itk::TimeProbe timer; // Timer for measuring the preprocessing time
	timer.Start();
	
	std::string write_file_path = file_path.substr(0, file_path.find_last_of('.'));
	
	try{
		Common::ReadImage3D(file_path, data_ptr);
		Common::WriteImage3D(write_file_path + std::string("_original.mhd"), data_ptr);
	}
	catch(itk::ExceptionObject& e){
		std::cout << e << std::endl;
		return EXIT_FAILURE;
	}

	//PARAMS: Make available outside function
	int median_radius = 1;
	int anis_diffusion_N_iter = 10;
	int anis_diffusion_conductance = 30;
	float smoothing_sigma = 2.0f;

	Common::MedianFilter(median_radius, data_ptr);
	Common::CurvatureAnisotropicDiffusion(anis_diffusion_N_iter, anis_diffusion_conductance, data_ptr);
	
	try{
		Common::WriteTIFFImage3D(write_file_path + std::string("_preprocessed.tif"), data_ptr);
	}
	catch(itk::ExceptionObject& e){
		std::cout << e << std::endl;
		return EXIT_FAILURE;
	}

	Common::GVFDiffusion(smoothing_sigma, write_file_path, data_ptr);	

	timer.Stop();
	std::cout << "The processing took " << timer.GetMeanTime() << " seconds. " << std::endl;

	return EXIT_SUCCESS;
}

void ftkVesselTracer::LoadPreprocessedData(std::string data_path){
	
	data_path = data_path.substr(0, data_path.find_last_of('.')); // get rid of the extension
	
	try{
		Common::ReadImage3D(data_path + std::string("_original.mhd"), this->originalData);
		Common::ReadImage3D(data_path + std::string(".mhd"), this->inputData); //input data is preprocessed data
		Common::ReadImage3D(data_path + std::string("_gx.mhd"), this->gx);
		Common::ReadImage3D(data_path + std::string("_gy.mhd"), this->gy);
		Common::ReadImage3D(data_path + std::string("_gz.mhd"), this->gz);
	}
	catch(itk::ExceptionObject& e){
		std::cout << e << std::endl;		
	}
}

void ftkVesselTracer::ComputeMaximumIntensityProjectionImage(void){
	
	typedef itk::MaximumProjectionImageFilter<RenderImageType3D, RenderImageType3D> MaxProjectionFilterType;
	MaxProjectionFilterType::Pointer max_intensity_projector = MaxProjectionFilterType::New();
	max_intensity_projector->SetInput(this->inputDataForRendering);

	max_intensity_projector->Update();

	this->MIPImage = max_intensity_projector->GetOutput();

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

void ftkVesselTracer::RenderMIPImage(void){

	ITKToVTKConnectorType::Pointer itk_to_vtk_connector = ITKToVTKConnectorType::New();
	
	//itk_to_vtk_connector->SetInput(max_intensity_image);
	itk_to_vtk_connector->SetInput(this->MIPImage);
	
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