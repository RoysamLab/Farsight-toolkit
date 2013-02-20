#include "Curvelet.h"


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
		_exit(0);
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






//int optionsCreate(const char* optfile, map<string,string>& options)
//{
//	options.clear();
//	ifstream fin(optfile); assert(fin.good());
//	string name;  fin>>name;
//	while(fin.good()) {
//		char cont[100];	 fin.getline(cont, 99);
//		options[name] = string(cont);
//		fin>>name;
//	}
//	fin.close();
//	return 0;
//}







int main(int argc, char** argv)
{
	clock_t ck0, ck1;
	ck0 = clock();
	time_t t1,t2;

	time(&t1);
	if(argc < 2 || argc > 3)
	{
		std::cout<<"Usage : curvelets.exe input_file sigma\n";
		//std::cout<<"Usage : curvelets.exe input_file '-options' [options_file]\n";
		return -1;
	}
	//assert(argc==3);
	//get options

	//map<string, string> opts;  
	//if(argc == 4)
	//	optionsCreate(argv[3], opts);

	//map<string,string>::iterator mi;

	//mi = opts.find("-nbangles_coarse"); 
	//if(mi!=opts.end())
	//{ istringstream ss((*mi).second); ss>>nbangles_coarse; }
	//else
	//{ nbangles_coarse = 8; printf("Chose nbangles_coarse = 8 as default\n");}

	//mi = opts.find("-nsigmas_coarse");
	//if(mi!=opts.end())
	//{ istringstream ss((*mi).second); ss>>nsigmas_coarse; }
	//else
	//{	  nsigmas_coarse = 2.2; printf("Chose nsigmas_coarse = 2.2 as default\n"); }

	//mi = opts.find("-nsigmas_fine"); 
	//if(mi!=opts.end())
	//{ istringstream ss((*mi).second); ss>>nsigmas_fine; }
	//else
	//{ nsigmas_fine = 2.5; printf("Chose nsigmas_fine = 2.5 as default\n");}

	//mi = opts.find("-neighb_weight"); 
	//if(mi!=opts.end())
	//{ istringstream ss((*mi).second); ss>>neighb_weight; }
	//else
	//{ neighb_weight = 0.5; printf("Chose neighb_weight = 0.5 as default\n"); }

	//mi = opts.find("-tuning_neighb"); 
	//if(mi!=opts.end())
	//{ istringstream ss((*mi).second); ss>>tuning_neighb; }
	//else
	//{ tuning_neighb = 0.6; printf("Chose tuning_neighb = 0.6 as default\n"); }
	//mi = opts.find("-sigma_ratio"); 
	//if (strcmp(argv[2], "-options"))
	//	{
	//		sigma_ratio = atof(argv[2]);
	//		std::cout << "sigma set to arg[2]";
	//	}
	//else if(mi!=opts.end())
	//{ istringstream ss((*mi).second); ss>>sigma_ratio; }
	//else
	//{ sigma_ratio = 0.2; printf("Chose sigma_ratio = 0.2 as default \n"); }
	//std::cout << "Sigma\t" << sigma_ratio << std::endl;
	//
	//mi = opts.find("-num_threads"); 
	//if(mi!=opts.end())
	//{ istringstream ss((*mi).second); ss>>numt; }
	//else
	//{ numt = 16; printf("Chose num_threads = 8 as default \n"); }

	//
	//mi = opts.find("-tile_size"); 
	//if(mi!=opts.end())
	//{ istringstream ss((*mi).second); ss>>tile_size; }
	//else
	//{ tile_size = 4024; printf("Chose tile_size = 1024 as default \n"); }

	//
	//
	//mi = opts.find("-border"); 
	//if(mi!=opts.end())
	//{ istringstream ss((*mi).second); ss>>border; }
	//else
	//{ border = 100; printf("Chose border = 100 as default \n"); }



	InputImageType::Pointer InputImage = readImage<InputImageType>(argv[1]);

	float sigma = atof(argv[2]);
	Curvelet curvletfilter = Curvelet();
	curvletfilter.SetSigma(sigma);
	InputImageType::Pointer outputim = curvletfilter.RunOnInputImage(InputImage);
	
	printf("writing the image to disk...\n");
	char buffer[1024];
	argv[1][strlen(argv[1])-4] = 0;
	sprintf(buffer, "%s_sigma_%4.2f_CV.tif",argv[1],sigma);
	writeImage<InputImageType>(outputim,buffer);	
	/*sprintf(buffer, "%s_CV_cos.mhd",argv[1]);
	writeImage<FloatImageType>(cosim,buffer);
	sprintf(buffer, "%s_CV_sin.mhd",argv[1]);
	writeImage<FloatImageType>(sinim,buffer);*/


	time(&t2);
	std::cout<<"Curvelet preprocessing takes "<< t2-t1 << " seconds (time_t calculation)\n";
	//ck1 = clock();  cout<<"Curvelet preprocessing takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds"<<endl;  ck0 = ck1;
	//scanf("%*d");
	return 0;
}

