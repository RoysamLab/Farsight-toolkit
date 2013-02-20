#include "helpers.h"
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
int main(int argc, char* argv[])
{
	if(argc <6)
	{
		std::cout<<"Usage: compute_extra_association <LabelImageFileName1> <LabelImageFileName2> <LabelImageFileName3> <ParametersFilename> <TableFilename>\n";
		return 0;
	}
	
	// read image files:
	OutputImageType::Pointer LabImageNuc = OutputImageType::New(); 
    LabImageNuc =  readImage<OutputImageType>(argv[1]);
	OutputImageType::Pointer LabImageDC; 
    LabImageDC =  readImage<OutputImageType>(argv[2]);
	OutputImageType::Pointer LabImageCA; 
    LabImageCA =  readImage<OutputImageType>(argv[3]);

	//read paramters file:
	std::map<std::string, std::string> opts; 
	 optionsCreate(argv[4], opts);
	 float xyr_nuc,zr_nuc,xyr_dc,zr_dc;
	//get input data
	 std::map<std::string,std::string>::iterator mi;
	mi = opts.find("xy_radius_nuc"); 
	if(mi!=opts.end()){ std::istringstream ss((*mi).second); ss>>xyr_nuc; }

	mi = opts.find("z_radius_nuc");
	if(mi!=opts.end()){ std::istringstream ss((*mi).second); ss>>zr_nuc; }
	
	mi = opts.find("xy_radius_dc");
	if(mi!=opts.end()){ std::istringstream ss((*mi).second); ss>>xyr_dc; }
	
	mi = opts.find("z_radius_dc");
	if(mi!=opts.end()){ std::istringstream ss((*mi).second); ss>>zr_dc; }

	printf("xy_radius_nuc:%f\n",xyr_nuc);
	printf("z_radius_nuc:%f\n",zr_nuc);
	printf("xy_radius_dc:%f\n",xyr_dc);
	printf("z_radius_dc:%f\n",zr_dc);

	//Load Table:
	vtkSmartPointer<vtkTable> table = NULL;
	if( file_exists(argv[5]) )
	{
		table = ftk::LoadTable(argv[5]);
		printf("loaded table...\n");
	}
	else
	{
		printf("could not read table...\n");
		return 0;
	}
	 //AnalyzeCAFeatures(LabImageNuc,LabImageCA,table);
	 AnalyzeNucDensity(LabImageNuc,xyr_nuc,zr_nuc,table);

	 ftk::SaveTable( argv[5], table );
	   

}