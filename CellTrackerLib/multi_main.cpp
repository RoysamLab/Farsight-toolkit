#include "helpers.h"
#include "MultiFrameCellTracker.h"
#include <ftkUtils.h>

using namespace helpers;
#if defined(_MSC_VER)
#pragma warning(disable: 4996)
#endif
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

	fvar.distVariance = 50;//8.77;//50
	fvar.distMean = 5;//3.3;//5
	fvar.spacing[0] = 1;
	fvar.spacing[1] = 1;
	fvar.spacing[2] = 4;
	fvar.timeVariance = 1;//.119//1;
	fvar.overlapVariance = 1;//0.034;//1;
	fvar.overlapMean = 0;//0.2;//0;
	fvar.variances[FeatureVariances::VOLUME] = 90000;//44000;//90000;
	fvar.MS_prior = 0.4;//
	fvar.AD_prior = 0.01;
	fvar.T_prior = 1;
	fvar.boundDistMean = 4;
	fvar.boundDistVariance = 50;

	if(argc==3)
	{
		 map<string, string> opts;  optionsCreate(argv[2], opts);
  
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
		//scanf("%*d");

	}
	else if(argc>=6)
	{
		printf("Usage: %s input_filenames labeled_filenames tracked_filenames [tracking_parameters_file]",argv[0]);
		scanf("%*d");
		return -1;
	}


	std::string line;
	// Read Input Filenames:
	std::vector< std::string > infnames;
	
	std::ifstream myfile (argv[1]);
	if (myfile.is_open())
	  {
		while ( myfile.good() )
		{
		  getline (myfile,line);
		  if(!line.empty())
		  {
			infnames.push_back(line);
		  }
		  cout << line << endl;
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
		  cout << line << endl;
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
		  cout << line << endl;
		}
		myfile2.close();
	  }
	else std::cout << "Unable to open tracked files"; 
	

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
	for(int t =0; t<num_time_points; t++)
	{
		tempimage = readImage<InputImageType>(infnames[t].c_str());	
		tempsegmented = readImage<LabelImageType>(labfnames[t].c_str());

		input_images.push_back(tempimage);
		label_images.push_back(tempsegmented);
	}

	// This will run the tracking :
	MultiFrameCellTracker * mfcellTracker = new MultiFrameCellTracker();
	mfcellTracker->set_parameters_from_cmd(fvar);
	mfcellTracker->set_inputs_from_cmd(input_images,label_images);
	std::vector< LabelImageType::Pointer > tracked_images = mfcellTracker->get_ouput_to_cmd();

	// Write the output:
	for(int t =0; t<tracked_images.size(); t++)
	{
		writeImage<LabelImageType>(tracked_images[t],trackfnames[t].c_str());
	}


	return 0;
}
