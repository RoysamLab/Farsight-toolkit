#include "helpers.h"
#include "MultiFrameCellTracker.h"
#include <ftkUtils.h>

using namespace helpers;
#if defined(_MSC_VER)
#pragma warning(disable: 4996)
#endif
template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
	printf("Writing %s ... \n",filename);
	typedef typename itk::ImageFileWriter<T> WriterType;

	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(im);
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

double start_t,end_t,diff_t;

bool file_exists(char *filename)
{
	FILE * fp = fopen(filename,"r");
	if(fp!=NULL)
	{
		fclose(fp);
		return true;
	}
	return false;
}


template <typename T>
typename T::Pointer readImage(const char *filename)
{
	printf("Reading %s ... \n",filename);
	typedef typename itk::ImageFileReader<T> ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();

	ReaderType::GlobalWarningDisplayOff();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		//return EXIT_FAILURE;
	}
	printf("Done\n");
	return reader->GetOutput();

}

int optionsCreate(const char* optfile, std::map<std::string,std::string>& options)
{
  options.clear();
  std::ifstream fin(optfile); assert(fin.good());
  std::string name;  fin>>name;
  while(fin.good()) {
	 char cont[100];	 fin.getline(cont, 99);
	 options[name] = std::string(cont);
	 fin>>name;
  }
  fin.close();
  return 0;
}






int main(int argc, char **argv)
{
	
	FeatureVariances fvar;
	for(int counter=0; counter < FeatureVariances::N; counter++)
	{
		fvar.variances[counter] = std::numeric_limits<float>::max();
	}

// 	fvar.distVariance = 50;//8.77;//50
// 	fvar.distMean = 5;//3.3;//5
// 	fvar.spacing[0] = 1;
// 	fvar.spacing[1] = 1;
// 	fvar.spacing[2] = 4;
// 	fvar.timeVariance = 1;//.119//1;
// 	fvar.overlapVariance = 1;//0.034;//1;
// 	fvar.overlapMean = 0;//0.2;//0;
// 	fvar.variances[FeatureVariances::VOLUME] = 90000;//44000;//90000;
// 	fvar.MS_prior = 0.4;//
// 	fvar.AD_prior = 0.01;
// 	fvar.T_prior = 1;
// 	fvar.boundDistMean = 4;
// 	fvar.boundDistVariance = 50;
// 	
// 	std::cout << std::endl << "HEREEEE " << fvar.FlagRoi << " ";

// 	if(argc==3)
// 	{
		 map<string, string> opts;  optionsCreate(argv[4], opts);
  
		//get input data
		map<string,string>::iterator mi;
		mi = opts.find("x_spacing"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.spacing[0]; }

		mi = opts.find("y_spacing");
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.spacing[1]; }
		
		mi = opts.find("z_spacing"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.spacing[2]; }
		
		mi = opts.find("distVariance");
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.distVariance; }

		mi = opts.find("distMean"); 
		if(mi!=opts.end())
		{ istringstream ss((*mi).second); ss>>fvar.distMean; }

		mi = opts.find("boundaryDistVariance");
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.boundDistVariance; }

		mi = opts.find("boundaryDistMean"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.boundDistMean; }

		mi = opts.find("timeVariance"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.timeVariance; }

		mi = opts.find("timeMean"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.timeMean; }

		mi = opts.find("overlapVariance"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.overlapVariance; }

		mi = opts.find("overlapMean"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.overlapMean; }

		mi = opts.find("volumeVariance");
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.variances[FeatureVariances::VOLUME]; }

		mi = opts.find("merge_split_prior"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.MS_prior; }

		mi = opts.find("appear_disappear_prior"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.AD_prior; }

		mi = opts.find("translation_prior"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.T_prior; }

		mi = opts.find("use_ROI"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.FlagRoi; }

		printf("fvar.spacing[0] = %f\n",fvar.spacing[0]);
		printf("fvar.spacing[1] = %f\n",fvar.spacing[1]);
		printf("fvar.spacing[2] = %f\n",fvar.spacing[2]);
		printf("fvar.distVariance = %f\n",fvar.distVariance);
		printf("fvar.distMean = %f\n",fvar.distMean);
		printf("fvar.boundDistVariance = %f\n",fvar.boundDistVariance);
		printf("fvar.boundDistMean = %f\n",fvar.boundDistMean);
		printf("fvar.timeVariance = %f\n", fvar.timeVariance);
		printf("fvar.timeMean = %f\n", fvar.timeMean);
		printf("fvar.overlapVariance = %f\n", fvar.overlapVariance);
		printf("fvar.overlapMean = %f\n", fvar.overlapMean);
		printf("fvar.volumeVariance = %f\n", fvar.variances[FeatureVariances::VOLUME]);
		printf("fvar.MS_prior = %f\n", fvar.MS_prior);
		printf("fvar.AD_prior = %f\n", fvar.AD_prior);
		printf("fvar.T_prior = %f\n", fvar.T_prior);
		printf("fvar.FlagRoi = %f\n", fvar.FlagRoi);
		//scanf("%*d");

// 	}
// 	else if(argc>=7)
// 	{
// 		printf("Usage: %s input_filenames labeled_filenames tracked_filenames [tracking_parameters_file]",argv[0]);
// 		scanf("%*d");
// 		return -1;
// 	}


	std::string line;
	// Read Input Filenames:
	std::vector< std::string > infnames;
	
	std::ifstream myfile (argv[1]);
	if (myfile.is_open())
	  {
		while ( myfile.good() )
		{
		  std::getline (myfile,line);
		  if(!line.empty())
		  {
			infnames.push_back(line);
		  }
		  std::cout << line << endl;
		}
		myfile.close();
	  }
	else std::cout << "Unable to open input files"; 
	// Read Label Filenames:
	std::vector< std::string > labfnames;
	
	std::ifstream myfile1 (argv[2]);
	if (myfile1.is_open())
	  {
		while ( myfile1.good() )
		{
		  getline (myfile1,line);
		  if(!line.empty())
		  {
			 labfnames.push_back(line);
		  }
		  std::cout << line << std::endl;
		}
		myfile1.close();
	  }
	else std::cout << "Unable to open labeled files"; 

	// Read Tracked Filenames:
	std::vector< std::string > trackfnames;
	
	std::ifstream myfile2 (argv[3]);
	if (myfile2.is_open())
	  {
		while ( myfile2.good() )
		{
		  getline (myfile2,line);
		  if(!line.empty())
		  {
			 trackfnames.push_back(line);
		  }
		  std::cout << line << std::endl;
		}
		myfile2.close();
	  }
	else std::cout << "Unable to open tracked files "<< argv[3] << " "; 
	
	// Read RoI Filenames:
	std::vector< std::vector< std::vector < int > > > wellsDetection;
	int numWells, numTimePoints;
	if( fvar.FlagRoi )
	{
		std::ifstream myfile3 (argv[5]);
// 		std::ifstream myfileOut ("/data/nicolas/test4/wellsRead.txt");
		if (myfile3.is_open())
		{
			while ( myfile3.good() )
			{
				myfile3 >> numWells >> numTimePoints;
				wellsDetection.resize(numWells);
				for(int i=0; i<numWells; ++i){
					wellsDetection.at(i).resize(numTimePoints);
					for(int j=0; j<numTimePoints; ++j){
						wellsDetection.at(i).at(j).resize(5);
						for(int k=0; k<5; ++k)
						{
							myfile3 >> wellsDetection.at(i).at(j).at(k);
						}
					}
				}
// 					std::cout << std::endl << "HOUston we have a problem: " << numWells <<" " << numTimePoints << " "<<argv[5];
// 				FILE* fp = fopen("/data/nicolas/test4/wellsRead.txt","w");
// 				for(int i=0; i<numWells; ++i){
// 					for(int j=0; j<numTimePoints; ++j){
// 						for(int k=0; k<5; ++k)
// 						{
// 							fprintf(fp,"%d\t",wellsDetection.at(i).at(j).at(k));
// 						}
// 					fprintf(fp,"\n");
// 					}
// 				}
// 				fclose(fp);			
			}
			myfile3.close();		
		}
		else std::cout << "Unable to open wells file files"<< argv[5] << " "; 
	}
	
	// Make sure this path is already created
	std::string resultsPath = argv[6];
	// All the other results which are not the final result
	std::string resultsPathTemp = argv[7];
	std::string resultsPathImage = argv[8];

	// Read Images:
	std::vector< InputImageType::Pointer > input_images;
	std::vector< LabelImageType::Pointer > label_images;

	int num_time_points = -1;
	if (infnames.size()==labfnames.size())
	{
		num_time_points = infnames.size();
	}
	else
	{
		std::cout<<infnames.size()<<std::endl;
		std::cout<<labfnames.size()<<std::endl;
		std::cout<< "The number of time points does not match (input,labels)\n";
		return -1;
	}

	LabelImageType::Pointer tempsegmented;
	InputImageType::Pointer tempimage;
	fvar.time_last = num_time_points;
// 	std::cout << std::endl << "IM HERE__" << num_time_points << " " << numWells << std::flush;
// 	scanf("%*d");
	
	for(int t =0; t<num_time_points; t++)
	{
		tempimage = readImage<InputImageType>(infnames[t].c_str());	
		tempsegmented = readImage<LabelImageType>(labfnames[t].c_str());

		input_images.push_back(tempimage);
		label_images.push_back(tempsegmented);
		
		stringstream ss;//create a stringstream
		ss << t << ".tif";
		string test = resultsPathImage+"/dTEST_Image_ " + ss.str();
		writeImage<InputImageType>(tempimage,test.c_str());
	}
	
	// All wells images
	std::vector< std::vector< LabelImageType::Pointer > > tracked_images;
	tracked_images.resize(numWells);
// 	std::vector< std::vector< LabelImageType::Pointer > > allWells;
// 	allWells.resize(numWells);
	
	// Min and Max Values
// 	std::vector< std::vector< unsigned short > > minAndMaxValues;
// 	minAndMaxValues.resize(numWells);
	
	// Maximum value of label of a particular well
	int maxValueofWell = 0;

	for(int i=0; i<numWells; ++i){

		// Cropping part
		std::vector< InputImageType::Pointer > input_images_well;
		std::vector< LabelImageType::Pointer > label_images_well;
		
		input_images_well.resize(num_time_points);
		label_images_well.resize(num_time_points);

		for(int j=0; j<num_time_points; ++j)
		{
			if( wellsDetection.at(i).at(j).at(3) == 1 ) // Square
			{
// 				std::cout << std::endl << "IM HERE__dude" << std::flush;
				
				int xCoord = wellsDetection.at(i).at(j).at(2);
				int yCoord = wellsDetection.at(i).at(j).at(1);
				
				int width = 32;
				int xMin = xCoord - width;
				int xMax = xCoord + width;
				int yMin = yCoord - width;
				int yMax = yCoord + width;

				// Create an image
				InputImageType::Pointer imageInputCropped = InputImageType::New();
				InputImageType::RegionType region;
				InputImageType::IndexType start;
				start[0] = 0;
				start[1] = 0;
				start[2] = 0;
				InputImageType::SizeType size;
				size[0] = 2*width+1;
				size[1] = 2*width+1;
				size[2] = 1;
				region.SetSize(size);
				region.SetIndex(start);
				imageInputCropped->SetRegions(region);
				imageInputCropped->Allocate();
				imageInputCropped->FillBuffer(0);
				try
				{
					imageInputCropped->Update();
				}
				catch(itk::ExceptionObject &err)
				{
					std::cerr << "ExceptionObject caught!" <<std::endl;
					std::cerr << err << std::endl;
				}
				
// 				std::cout << std::endl << "IM HERE__dude2" << std::flush;
				
				LabelImageType::Pointer imageLabelCropped = LabelImageType::New();
				LabelImageType::RegionType region2;
				LabelImageType::IndexType start2;
				start2[0] = 0;
				start2[1] = 0;
				start2[2] = 0;
				LabelImageType::SizeType size2;
				size2[0] = 2*width+1;
				size2[1] = 2*width+1;
				size2[2] = 1;
				region2.SetSize(size2);
				region2.SetIndex(start2);
				imageLabelCropped->SetRegions(region2);
				imageLabelCropped->Allocate();
				imageLabelCropped->FillBuffer(0);
				try
				{
					imageLabelCropped->Update();
				}
				catch(itk::ExceptionObject &err)
				{
					std::cerr << "ExceptionObject caught!" <<std::endl;
					std::cerr << err << std::endl;
				}
				
// 				std::cout << std::endl << "IM HERE__dude3" << std::flush;
				InputImageType::PixelType * inputImagePointer = input_images.at(j)->GetBufferPointer();
				LabelImageType::PixelType * labelImagePointer = label_images.at(j)->GetBufferPointer();
				
				InputImageType::PixelType * inputCroppedImagePointer = imageInputCropped->GetBufferPointer();
				LabelImageType::PixelType * labelCroppedImagePointer = imageLabelCropped->GetBufferPointer();
				
				for( int xx = xMin; xx <= xMax; ++xx )
				{
					for( int yy = yMin; yy <= yMax; ++yy )
					{
						//Boundary conditions
						int indexOriginal = xx+yy*512;
						int indexCropped = (xx-xMin)+(yy-yMin)*(2*width+1);
						if( (xx<0) || (yy<0) || (xx>=512) || (yy>=512) )
						{
							inputCroppedImagePointer[indexCropped] = 0;
							labelCroppedImagePointer[indexCropped] = 0;
						}
						else
						{
							inputCroppedImagePointer[indexCropped] = inputImagePointer[indexOriginal];
							labelCroppedImagePointer[indexCropped] = labelImagePointer[indexOriginal];
						}
					}
				}
// 				std::cout << std::endl << "IM HERE__dude4" << num_time_points << std::flush;
				input_images_well.at(j) = imageInputCropped;
// 				std::cout << std::endl << "IM HERE__dude6" << std::flush;
				label_images_well.at(j) = imageLabelCropped;
// 				std::cout << std::endl << "IM HERE__dude5" << std::flush;
			}
			else if( wellsDetection.at(i).at(j).at(3) == 2 ) // Square
			{
// 				std::cout << std::endl << "IM HERE_romb__dude" << std::flush;	
				
				int xCoord = wellsDetection.at(i).at(j).at(1);
				int yCoord = wellsDetection.at(i).at(j).at(2);
				
				int width = 42;
				int xMin = xCoord - width;
				int xMax = xCoord + width;
				int yMin = yCoord - width;
				int yMax = yCoord + width;

				// Create an image
				InputImageType::Pointer imageInputCropped = InputImageType::New();
				InputImageType::RegionType region;
				InputImageType::IndexType start;
				start[0] = 0;
				start[1] = 0;
				start[2] = 0;
				InputImageType::SizeType size;
				size[0] = 2*width+1;
				size[1] = 2*width+1;
				size[2] = 1;
				region.SetSize(size);
				region.SetIndex(start);
				imageInputCropped->SetRegions(region);
				imageInputCropped->Allocate();
				imageInputCropped->FillBuffer(0);
				try
				{
					imageInputCropped->Update();
				}
				catch(itk::ExceptionObject &err)
				{
					std::cerr << "ExceptionObject caught!" <<std::endl;
					std::cerr << err << std::endl;
				}
// 				imageInputCropped->FillBuffer(0);
				
// 				std::cout << std::endl << "IM HERE__dude2" << std::flush;
				
				LabelImageType::Pointer imageLabelCropped = LabelImageType::New();
				LabelImageType::RegionType region2;
				LabelImageType::IndexType start2;
				start2[0] = 0;
				start2[1] = 0;
				start2[2] = 0;
				LabelImageType::SizeType size2;
				size2[0] = 2*width+1;
				size2[1] = 2*width+1;
				size2[2] = 1;
				region2.SetSize(size2);
				region2.SetIndex(start2);
				imageLabelCropped->SetRegions(region2);
				imageLabelCropped->Allocate();
				imageLabelCropped->FillBuffer(0);
				try
				{
					imageLabelCropped->Update();
				}
				catch(itk::ExceptionObject &err)
				{
					std::cerr << "ExceptionObject caught!" <<std::endl;
					std::cerr << err << std::endl;
				}
				
				InputImageType::PixelType * inputImagePointer = input_images.at(j)->GetBufferPointer();
				LabelImageType::PixelType * labelImagePointer = label_images.at(j)->GetBufferPointer();
				
				InputImageType::PixelType * inputCroppedImagePointer = imageInputCropped->GetBufferPointer();
				LabelImageType::PixelType * labelCroppedImagePointer = imageLabelCropped->GetBufferPointer();
				
				for( int xx = xMin; xx <= xMax; ++xx )
				{
					for( int yy = yMin; yy <= yMax; ++yy )
					{
						//Boundary conditions
						int indexOriginal = xx+yy*512;
						int indexCropped = (xx-xMin)+(yy-yMin)*(2*width+1);
						if( (xx<0) || (yy<0) || (xx>=512) || (yy>=512) )
						{
							inputCroppedImagePointer[indexCropped] = 0;
							labelCroppedImagePointer[indexCropped] = 0;
						}
						else if( std::abs(xx-xMin) + std::abs(yy-yMin) < width + 1 )
						{
							inputCroppedImagePointer[indexCropped] = inputImagePointer[indexOriginal];
							labelCroppedImagePointer[indexCropped] = labelImagePointer[indexOriginal];
						}
						else 
						{
							inputCroppedImagePointer[indexCropped] = 0;
							labelCroppedImagePointer[indexCropped] = 0;
						}
					}
				}
				input_images_well.at(j) = imageInputCropped;
				label_images_well.at(j) = imageLabelCropped;
			}
// 			stringstream ss;//create a stringstream
// 			ss << i << "__" << j << ".tif";
// 			string test = resultsPathTemp+"/cTEST_Input_ " + ss.str();
// 			writeImage<InputImageType>(input_images_well.at(j),test.c_str());
// 			
// 			stringstream ss2;//create a stringstream
// 			ss2 << i << "__" << j << ".tif";
// 			string test2 = resultsPathTemp+"/cTEST_Label_ " + ss2.str();
// 			writeImage<LabelImageType>(label_images_well.at(j),test2.c_str());
		}

		// This will run the tracking :
		MultiFrameCellTracker * mfcellTracker = new MultiFrameCellTracker();
		mfcellTracker->set_parameters_from_cmd(fvar);
		//
		mfcellTracker->set_inputs_from_cmd(input_images_well,label_images_well);
		tracked_images.at(i) = mfcellTracker->get_ouput_to_cmd();
		

// 		for(int j=0; j<num_time_points; ++j)
// 		{
// 		
// 			stringstream ss2;//create a stringstream
// 			ss2 << i << "__" << j << ".tif";
// 			string test2 = resultsPathTemp+"/bTrackResult_ " + ss2.str();
// 			writeImage<LabelImageType>(tracked_images.at(i).at(j),test2.c_str());
// 			
// 		}
		
		maxValueofWell = relabelWells(tracked_images.at(i),maxValueofWell);
		
		
		delete mfcellTracker;
	}
	
	std::vector< LabelImageType::Pointer > allBigLabelImages;
	allBigLabelImages.resize(num_time_points);
	
	
// 	std::cout << std::endl << "IM HERE__dude__9909" << std::flush;
	
	for(int j=0; j<num_time_points; ++j)
	{
		LabelImageType::Pointer imageLabelBig = LabelImageType::New();
		LabelImageType::RegionType region2;
		LabelImageType::IndexType start2;
		start2[0] = 0;
		start2[1] = 0;
		start2[2] = 0;
		LabelImageType::SizeType size2;
		size2[0] = 512;
		size2[1] = 512;
		size2[2] = 1;
		region2.SetSize(size2);
		region2.SetIndex(start2);
		imageLabelBig->SetRegions(region2);
		imageLabelBig->Allocate();
		imageLabelBig->FillBuffer(0);
		try
		{
			imageLabelBig->Update();
		}
		catch(itk::ExceptionObject &err)
		{
			std::cerr << "ExceptionObject caught!" <<std::endl;
			std::cerr << err << std::endl;
		}
		LabelImageType::PixelType * labelImageBigPointer = imageLabelBig->GetBufferPointer();
		
		for( int i=0;i<numWells;++i )
		{
			if( wellsDetection.at(i).at(j).at(3) == 1 ) // Square
			{
// 				std::cout << std::endl << "IM HERE__dude" << std::flush;
				
				int xCoord = wellsDetection.at(i).at(j).at(2);
				int yCoord = wellsDetection.at(i).at(j).at(1);
				
				int width = 32;
				int xMin = xCoord - width;
				int xMax = xCoord + width;
				int yMin = yCoord - width;
				int yMax = yCoord + width;
				
				LabelImageType::PixelType * labelImageWellPointer = tracked_images.at(i).at(j)->GetBufferPointer();

				for( int xx = xMin; xx <= xMax; ++xx )
				{
					for( int yy = yMin; yy <= yMax; ++yy )
					{
						//Boundary conditions
						int indexOriginal = xx+yy*512;
						int indexCropped = (xx-xMin)+(yy-yMin)*(2*width+1);
						if( !( (xx<0) || (yy<0) || (xx>=512) || (yy>=512) ) )
						{
							labelImageBigPointer[indexOriginal] = labelImageWellPointer[indexCropped];
						}
					}
				}
			}

					
			else if( wellsDetection.at(i).at(j).at(3) == 2 ) // Square
			{
// 				std::cout << std::endl << "IM HERE_romb__dude" << std::flush;	
				
				int xCoord = wellsDetection.at(i).at(j).at(1);
				int yCoord = wellsDetection.at(i).at(j).at(2);
				
				int width = 42;
				int xMin = xCoord - width;
				int xMax = xCoord + width;
				int yMin = yCoord - width;
				int yMax = yCoord + width;
				
				LabelImageType::PixelType * labelImageWellPointer = tracked_images.at(i).at(j)->GetBufferPointer();

				for( int xx = xMin; xx <= xMax; ++xx )
				{
					for( int yy = yMin; yy <= yMax; ++yy )
					{
						//Boundary conditions
						int indexOriginal = xx+yy*512;
						int indexCropped = (xx-xMin)+(yy-yMin)*(2*width+1);
						if( !( (xx<0) || (yy<0) || (xx>=512) || (yy>=512) ) )
						{
							// Dont use the condition of the rombus since the tracking is not supposed to change the position of the cells, only the labels 
							// else if( std::abs(xx-xMin) + std::abs(yy-yMin) < width + 1 )
							labelImageBigPointer[indexOriginal] = labelImageWellPointer[indexCropped];
						}
					}
				}
			}
		}

		std::stringstream ss3;//create a stringstream
		ss3 << j << "__" << ".tif";
		std::string test3 = resultsPath+"/aRelabelTrackResult_" + ss3.str();
		writeImage<LabelImageType>(imageLabelBig,test3.c_str());
	}
	
	
	


// 	// Read Images:
// 	std::vector< InputImageType::Pointer > input_images;
// 	std::vector< LabelImageType::Pointer > label_images;
// 
// 	int num_time_points = -1;
// 	if (infnames.size()==labfnames.size())
// 	{
// 		num_time_points = infnames.size();
// 	}
// 	else
// 	{
// 		std::cout<<infnames.size()<<std::endl;
// 		std::cout<<labfnames.size()<<std::endl;
// 		std::cout<< "The number of time points does not match (input,labels)\n";
// 		return -1;
// 	}
// 
// 	LabelImageType::Pointer tempsegmented;
// 	InputImageType::Pointer tempimage;
// 	fvar.time_last = num_time_points;
// 	for(int t =0; t<num_time_points; t++)
// 	{
// 		tempimage = readImage<InputImageType>(infnames[t].c_str());	
// 		tempsegmented = readImage<LabelImageType>(labfnames[t].c_str());
// 
// 		input_images.push_back(tempimage);
// 		label_images.push_back(tempsegmented);
// 	}
// 
// 	// This will run the tracking :
// 	MultiFrameCellTracker * mfcellTracker = new MultiFrameCellTracker();
// 	mfcellTracker->set_parameters_from_cmd(fvar);
// 	//
// 	mfcellTracker->set_inputs_from_cmd(input_images,label_images);
// 	std::vector< LabelImageType::Pointer > tracked_images = mfcellTracker->get_ouput_to_cmd();
// 
// 	// Write the output:
// 	for(int t =0; t<tracked_images.size(); t++)
// 	{
// 		writeImage<LabelImageType>(tracked_images[t],trackfnames[t].c_str());
// 	}


	return 0;
}
