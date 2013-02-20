//includes

#include <stdio.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <iomanip>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkMedianImageFilter.h>
#include <itkMeanImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkVector.h>
#include <itkVTKImageExport.h>
#include <itkVTKImageImport.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_hungarian_algorithm.h>
#include <vnl/algo/vnl_qr.h>
//#include <ilcplex/ilocplex.h>
//ILOSTLBEGIN

#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))

//typedefs

typedef unsigned char InputPixelType;
typedef unsigned char OutputPixelType;

typedef itk::Vector<unsigned char, 3> VectorPixelType;
typedef itk::Image<VectorPixelType, 3> ColorImageType;
typedef itk::Image<VectorPixelType, 2> Color2DImageType;

typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<OutputPixelType,3> OutputImageType;
typedef itk::Image<short int,3> LabelImageType;

typedef itk::Image<InputPixelType,2> Input2DImageType;
typedef itk::Image<InputPixelType,2> Output2DImageType;

typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
typedef itk::ImageRegionIterator<InputImageType> IteratorType;

typedef itk::ImageRegionConstIterator<LabelImageType> ConstLabelIteratorType;
typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;


typedef itk::ImageRegionConstIterator<Input2DImageType> Const2DIteratorType;
typedef itk::ImageRegionIterator<Input2DImageType> twoDIteratorType;

typedef itk::ImageRegionConstIterator<ColorImageType> ConstColorIteratorType;
typedef itk::ImageRegionIterator<ColorImageType> ColorIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Input2DImageType > LinearIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType > SliceIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Color2DImageType > LinearColorIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ColorImageType > SliceColorIteratorType;

typedef itk::ConstNeighborhoodIterator<InputImageType> NeighborhoodIteratorType;

typedef itk::MedianImageFilter<InputImageType,InputImageType> MedianFilterType;
typedef itk::BinaryMedianImageFilter<InputImageType,InputImageType> BinaryMedianFilterType;

typedef itk::Image<bool,3> BoolImageType;
typedef itk::BinaryThresholdImageFilter<InputImageType,OutputImageType> ThresholdFilterType;
typedef itk::OtsuThresholdImageFilter<InputImageType,OutputImageType> OtsuThresholdFilterType;

typedef itk::Image<short int,3> DistanceImageType;
typedef itk::DanielssonDistanceMapImageFilter<InputImageType,DistanceImageType> DistanceMapFilterType;
typedef DistanceMapFilterType::VectorImageType OffsetImageType;

typedef itk::ConnectedComponentImageFilter<InputImageType,LabelImageType> ConnectedFilterType;
typedef itk::RelabelComponentImageFilter<LabelImageType,LabelImageType> RelabelFilterType;




//structures

struct Feature{
	int x,y,z,t;
	int volume,tag;
	int num;
	double diff(Feature &other)
	{
		double sum;
		sum = 1-exp(-double((x-other.x)*(x-other.x)+(y-other.y)*(y-other.y)+(z-other.z)*(z-other.z))/2/25);
		return sum;
	}
};

//functions



/**
* This function will connect the given itk::VTKImageExport filter to
* the given vtkImageImport filter.
*/
template <typename ITK_Exporter, typename VTK_Importer>
void ConnectPipelines(ITK_Exporter exporter, VTK_Importer* importer)
{
	importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
	importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
	importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
	importer->SetSpacingCallback(exporter->GetSpacingCallback());
	importer->SetOriginCallback(exporter->GetOriginCallback());
	importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
	importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
	importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
	importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
	importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
	importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
	importer->SetCallbackUserData(exporter->GetCallbackUserData());
}


template <typename T>
typename T::Pointer readImage(const char *filename)
{
	printf("Reading %s ... ",filename);
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
		std::cout << "ExceptionObject caught!" <<std::endl;
		std::cout << err << std::endl;
		//return EXIT_FAILURE;
	}
	printf("Done.\n");
	return reader->GetOutput();

}
template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
	printf("Writing %s ... ",filename);
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
		std::cout << "ExceptionObject caught!" <<std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}
	printf("Done.\n");
	return EXIT_SUCCESS;
}

template <typename T>
typename T::Pointer getEmpty(typename T::Pointer input)
{
	typename T::Pointer output = T::New();
	output->SetRegions(input->GetLargestPossibleRegion());
	output->Allocate();
	output->FillBuffer(0);
	output->Update();
	return output;
}

InputImageType::Pointer getEmpty(int s1,int s2, int s3)
{
	InputImageType::Pointer p = InputImageType::New();
	InputImageType::SizeType size;
	InputImageType::IndexType index;
	InputImageType::RegionType region;
	size[0] = s1; size[1] = s2; size[2] = s3;
	index.Fill(0);
	region.SetSize(size);
	region.SetIndex(index);
	p->SetRegions(region);
	p->Allocate();
	return p;
}



template <typename T3D, typename T2D>
typename T2D::Pointer getProjection(typename T3D::Pointer im)
{
	typename T2D::Pointer output = T2D::New();
	typename T2D::RegionType region;
	typename T2D::SizeType size;
	typename T2D::IndexType index;
	index[0]=0;
	index[1]=0;
	size[0] = im->GetLargestPossibleRegion().GetSize()[0];
	size[1] = im->GetLargestPossibleRegion().GetSize()[1];
	region.SetSize(size);
	region.SetIndex(index);
	output->SetRegions(region);
	output->Allocate();

	SliceIteratorType inputIt(im,im->GetLargestPossibleRegion());
	LinearIteratorType outputIt(output,output->GetLargestPossibleRegion());

	inputIt.SetFirstDirection(0);
	inputIt.SetSecondDirection(1);
	outputIt.SetDirection(0);


	outputIt.GoToBegin();
	while ( ! outputIt.IsAtEnd() )
	{
		while ( ! outputIt.IsAtEndOfLine() )
		{
			outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
			++outputIt;
		}
		outputIt.NextLine();
	}

	inputIt.GoToBegin();
	outputIt.GoToBegin();

	while( !inputIt.IsAtEnd() )
	{
		while ( !inputIt.IsAtEndOfSlice() )
		{
			while ( !inputIt.IsAtEndOfLine() )
			{
				outputIt.Set( MAX( outputIt.Get(), inputIt.Get() ));
				++inputIt;
				++outputIt;
			}
			outputIt.NextLine();
			inputIt.NextLine();
		}
		outputIt.GoToBegin();
		inputIt.NextSlice();

	}
	return output;
}

Input2DImageType::Pointer getCollage(InputImageType::Pointer im[4])
{
	InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
	Input2DImageType::Pointer imcollage = Input2DImageType::New();
	Input2DImageType::SizeType size2d;
	Input2DImageType::RegionType region;
	Input2DImageType::IndexType index;
	size2d[0] = size[0]*2;
	size2d[1] = size[1]*2;
	index[0]=0;
	index[1]=0;

	region.SetIndex(index);
	region.SetSize(size2d);
	imcollage->SetRegions(region);
	imcollage->Allocate(); 

	int startpoints[][2]={0,0,size[0],0,0,size[1],size[0],size[1]};
	size2d[0]=size[0];
	size2d[1]=size[1];
	region.SetSize(size2d);
	for(int co = 0; co<4; co++)
	{
		index[0]=startpoints[co][0];
		index[1]=startpoints[co][1];
		region.SetIndex(index);
		Input2DImageType::Pointer proj = getProjection<InputImageType, Input2DImageType>(im[co]);
		Const2DIteratorType inputIt(proj,proj->GetLargestPossibleRegion());
		twoDIteratorType outputIt(imcollage,region);
		for(inputIt.GoToBegin(),outputIt.GoToBegin();!inputIt.IsAtEnd();++inputIt,++outputIt)
		{
			outputIt.Set(inputIt.Get());
		}
	}

	return imcollage;

}


Input2DImageType::Pointer getGalleryView(InputImageType::Pointer im[], int rows, int cols)
{
	InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
	Input2DImageType::Pointer imcollage = Input2DImageType::New();
	Input2DImageType::SizeType size2d;
	Input2DImageType::RegionType region;
	Input2DImageType::IndexType index;
	size2d[0] = size[0]*cols;
	size2d[1] = size[1]*rows;
	index[0]=0;
	index[1]=0;

	region.SetIndex(index);
	region.SetSize(size2d);
	imcollage->SetRegions(region);
	imcollage->Allocate(); 

	int startpoints[2000][2];
	for(int co = 0; co < rows*cols ; co++)
	{
		startpoints[co][0] = size[0]*(co%cols);
		startpoints[co][1] = size[1]*(co/cols);
	}

	size2d[0]=size[0];
	size2d[1]=size[1];
	region.SetSize(size2d);
	for(int co = 0; co<rows*cols; co++)
	{
		index[0]=startpoints[co][0];
		index[1]=startpoints[co][1];
		region.SetIndex(index);
		Input2DImageType::Pointer proj = getProjection<InputImageType, Input2DImageType>(im[co]);
		Const2DIteratorType inputIt(proj,proj->GetLargestPossibleRegion());
		twoDIteratorType outputIt(imcollage,region);
		for(inputIt.GoToBegin(),outputIt.GoToBegin();!inputIt.IsAtEnd();++inputIt,++outputIt)
		{
			outputIt.Set(inputIt.Get());
		}
	}

	return imcollage;
}


vnl_matrix<double> getFingerPrintMatrix()
{
	FILE *fp = fopen("L:\\unmixing\\SinglePositives\\spectral\\estimated_matrix.txt","r");
	int nrows = 15;
	int ncols = 6;
	
	//int arr[] = {0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0};
	int arr[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	//int arr[] = {0, 1, 0, 1, 1, 1};
	int sum = 0;
	for(int co = 0; co < nrows; co++)
	{
		sum = sum + arr[co];
	}
	vnl_matrix<double> start(sum,ncols);
	for(int coy = 0; coy < ncols; coy ++)
	{
		int pc = 0;
		for(int cox =0; cox < nrows; cox ++)
		{
			double temp;
			fscanf(fp,"%lf",&temp);
			if(arr[cox] == 1)
			{
				start[pc++][coy] = temp;
			}
		}
	}
	fclose(fp);
	return start;
}
vnl_matrix<double> getFingerPrintMatrix_defunct()
{
#define BASE_FOLDER  "C:\\Users\\arun\\Research\\unmixing_work\\data\\"
	
	vnl_matrix<double> mat(16,4);
	float sum = 0;
	FILE *fp = fopen(BASE_FOLDER "2Map2TauSinglePositive.txt","r");
	
	for(int counter =0; counter < 15; counter++)
	{
		int t1; float t2;
		fscanf(fp, "%d %f",&t1, &t2);
		mat(counter+1,0)= t2;
		sum = sum + t2*t2;
	}
	mat(0,0) = mat(1,0);
	sum = sum + mat(1,0)*mat(1,0);
	sum = sqrt(sum);
	for(int counter = 0; counter < 16; counter++)
		mat(counter,0) = mat(counter,0)/sum;
	fclose(fp);

	printf("Done one file ....\n");
	sum=0;
	fp = fopen(BASE_FOLDER "CNPase1300.txt","r");
	
	for(int counter =0; counter < 15; counter++)
	{
		int t1; float t2;
		fscanf(fp, "%d %f",&t1, &t2);
		mat(counter+1,1)= t2;
		sum = sum + t2*t2;
	}
	mat(0,1) = mat(1,1);
	sum = sum + mat(1,1)*mat(1,1);
	sum = sqrt(sum);
	for(int counter = 0; counter < 16; counter++)
		mat(counter,1) = mat(counter,1)/sum;
	fclose(fp);
printf("Done one file ....\n");
	sum=0;
	fp = fopen(BASE_FOLDER "CyQuantSinglePositive.txt","r");
	
	for(int counter =0; counter < 15; counter++)
	{
		int t1; float t2;
		fscanf(fp, "%d %f",&t1, &t2);
		mat(counter+1,2)= t2;
		sum = sum + t2*t2;
	}
	mat(0,2) = mat(1,2);
	sum = sum + mat(1,2)*mat(1,2);
	sum = sqrt(sum);
	for(int counter = 0; counter < 16; counter++)
		mat(counter,2) = mat(counter,2)/sum;
	fclose(fp);
printf("Done one file ....\n");
	sum=0;
	fp = fopen(BASE_FOLDER "LectinSinglePositive2.txt","r");
	
	for(int counter =0; counter < 15; counter++)
	{
		int t1; float t2;
		fscanf(fp, "%d %f",&t1, &t2);
		mat(counter+1,3)= t2;
		sum = sum + t2*t2;
	}
	mat(0,3) = mat(1,3);
	sum = sum + mat(1,3)*mat(1,3);
	sum = sqrt(sum);
	for(int counter = 0; counter < 16; counter++)
		mat(counter,3) = mat(counter,3)/sum;
	fclose(fp);
printf("Done one file ....\n");
	return mat;
	
}


double square_function (double a)
{
	return a*a;
}
void unmix_cluster_from_matrix(vnl_matrix<double> mixed, vnl_matrix<double> &start, vnl_vector<unsigned char> &indices)
{
	printf("Entered unmix_clsuter_from_matrix\n");
	int m = start.cols();
	int n = start.rows();

	bool converge = false;

	start.normalize_columns();
	// some neat way of generating sum(k.^2,2)*ones(1,n);
	vnl_matrix<double> sqmixed(mixed);
	sqmixed.apply(square_function);
	vnl_matrix<double> ones(n,1);
	ones.fill(1);
	sqmixed = sqmixed*ones;
	ones = ones.transpose();
	sqmixed = sqmixed*ones;
	// end of neat way
	vnl_vector<double> max_along_rows(mixed.rows());


	int iteration_number = 0;
	vnl_matrix<double> old_start;
	start.print(std::cout);
	printf("Entered the loop\n");
	while(!converge)
	{
		iteration_number++;
		if(iteration_number>50)
			break;
		vnl_matrix<double> out = mixed*start;
		//printf("out.size = %d %d\n",out.rows(),out.columns());
	//	out.apply(square_function);
	//	out -= sqmixed;
		max_along_rows.fill(-16466445);
		int num_rows = mixed.rows();
		int num_columns = m;
		old_start = start;
		start.fill(0);// set the current estimate to a zero matrix.
		//printf("Finding maximum along components\n");
		printf("-%d\n", iteration_number);
		for(int counterx = 0; counterx < num_rows; counterx++)
		{
			//printf("Before first loop\n");
			indices[counterx]=-1;
			for(int countery = 0; countery < num_columns; countery++)
			{
				if(max_along_rows[counterx] < (out[counterx][countery]))
				{
					max_along_rows[counterx] = (out[counterx][countery]);
					indices[counterx] = countery;
				}
			}
			//printf("After first loop\n");
			double temp=0;
			for(int countery = 0; countery < n; countery++)
			{
			//	printf("start[%d][%d] += out[%d][%d] start.size %d %d\n",countery,indices[counterx],counterx,countery,start.rows(),start.columns());
			//	printf("start->data = %p\n",start.data_block());
				start[countery][indices[counterx]] += mixed[counterx][countery];
				
				
			}
			//printf("After Second loop\n");
		}
		start.normalize_columns();
		start.print(std::cout);

		//printf("About to find the change\n");
	 	double change = sqrt((old_start-start).apply(square_function).absolute_value_sum());
		if(change < 0.0001)
			converge = true;
		//printf("iteration %d complete\n",iteration_number);
	}
}

void unmix_cluster_from_matrix_v2(vnl_matrix<double> mixed, vnl_matrix<double> &start, vnl_vector<unsigned char> &indices)
{
	printf("Entered unmix_clsuter_from_matrix\n");
	int m = start.cols();
	int n = start.rows();

	std::ofstream f;
	f.open("C:/Users/arun/Research/unmixing_work/build/release/matrix.txt",std::ios_base::out);
	bool converge = false;

	start.normalize_columns();
	// some neat way of generating sum(k.^2,2)*ones(1,n);
	vnl_matrix<double> sqmixed(mixed);
	sqmixed.apply(square_function);
	vnl_matrix<double> ones(n,1);
	ones.fill(1);
	sqmixed = sqmixed*ones;
	ones = ones.transpose();
	sqmixed = sqmixed*ones;
	// end of neat way
	vnl_vector<double> max_along_rows(mixed.rows());

	std::vector< vnl_vector<double> > subsets[10];
	int iteration_number = 0;
	vnl_matrix<double> old_start;
	vnl_matrix<double> A[10];
	start.print(std::cout);
	printf("Entered the loop\n");
	while(!converge)
	{
		iteration_number++;
		if(iteration_number>50)
			break;
		vnl_matrix<double> out = mixed*start;
		//printf("out.size = %d %d\n",out.rows(),out.columns());
	//	out.apply(square_function);
	//	out -= sqmixed;
		max_along_rows.fill(-16466445);
		int num_rows = mixed.rows();
		int num_columns = m;
		old_start = start;
		start.fill(0);// set the current estimate to a zero matrix.
		for(int counter=0; counter< m; counter++)
		{
			subsets[counter].clear();
		}
		//printf("Finding maximum along components\n");
		printf("-%d\n", iteration_number);
		for(int counterx = 0; counterx < num_rows; counterx++)
		{
			//printf("Before first loop\n");
			indices[counterx]=-1;
			for(int countery = 0; countery < num_columns; countery++)
			{
				if(max_along_rows[counterx] < out[counterx][countery])
				{
					max_along_rows[counterx] = out[counterx][countery];
					indices[counterx] = countery;
				}
			}
			//printf("After first loop\n");
			double temp=0;

			subsets[indices[counterx]].push_back(mixed.get_row(counterx));
			for(int countery = 0; countery < n; countery++)
			{
			//	printf("start[%d][%d] += out[%d][%d] start.size %d %d\n",countery,indices[counterx],counterx,countery,start.rows(),start.columns());
			//	printf("start->data = %p\n",start.data_block());
				//start[countery][indices[counterx]] += mixed[counterx][countery];
				
				
				
			}
			//printf("After Second loop\n");
		}
		
		for(int counter = 0; counter < m; counter++)
		{
			vnl_vector<double> p(n);
			p.fill(1);
			vnl_vector<double> t(n),old_p(n);
			t.fill(0);
			old_p = p;
			while(1)
			{
				t.fill(0);
				for(int counter1 = 0; counter1 < subsets[counter].size(); counter1++)
				{
					t = t + (dot_product(subsets[counter][counter1],p)*subsets[counter][counter1]);
				}
				p = t.normalize();
				if((p-old_p).two_norm()< 1e-6)
					break;
				old_p = p;
			}
			start.set_column(counter,p);
		}
		start.normalize_columns();
		start.print(std::cout);
		if(start.get_column(0).is_zero() == 1 || start.get_column(1).is_zero() == 1)
		{
			srand(time(NULL));
			printf("resetting...\n");
			for(int counter=0; counter<m; counter++)
			{
				for(int ncounter = 0; ncounter < n; ncounter++)
				{
					start[ncounter][counter] = rand()%100;
				}
			}
			std::cout<<std::setprecision(2);
			start.print(std::cout);
			start.normalize_columns();
			continue;
		}

		//printf("About to find the change\n");
	 	double change = (old_start-start).apply(square_function).absolute_value_sum();
		if(change < 0.00001)
			converge = true;
		//printf("iteration %d complete\n",iteration_number);
	}
	start.print(f);
	f.close();

}

#define DO_FILTERING 1
//#define NO_SATURATED_PIXELS_IN_UNMIXIMG
#define MAX_CHANNELS 50
#define NUM_ELEMENTS_FOR_UNMIXING 1000000
#define MIN_NORM 30
#define UNMIX_CLUSTERING_PIXEL_CLASSIFICATION
//#define UNMIX_METHOD3 
//#define NORMALIZE_OUTPUT_ROWS 
//#define UNMIX_METHOD4
//#define GENERATE_FINGERPRINT_MODE

void my_normalize_rows(vnl_matrix<double> &mat)
{
	for(int co = 0; co < mat.rows(); co++)
	{
		vnl_vector<double> r = mat.get_row(co);
		double sum = r.one_norm();
		r = r/sum;
		mat.set_row(co,r);
	}
}

void my_normalize_columns(vnl_matrix<double> &mat)
{
	for(int co = 0; co < mat.columns(); co++)
	{
		vnl_vector<double> r = mat.get_column(co);
		double sum = r.one_norm();
		r = r/sum;
		mat.set_column(co,r);
	}
}

void unmix_clustering(InputImageType::Pointer im[],InputImageType::Pointer om[],int n, int m)
{

	
	printf("Unnmixing %d channels ...\n",n);
	InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();

	printf("I'm here\n");
	MedianFilterType::Pointer filt[MAX_CHANNELS];
	typedef itk::MeanImageFilter<InputImageType,InputImageType> MeanFilterType;
	MeanFilterType::Pointer filt1[MAX_CHANNELS];
	IteratorType iterator[MAX_CHANNELS];
	IteratorType assigniter[MAX_CHANNELS];

	InputImageType::SizeType radius;


	InputImageType::SizeType imagesize = im[0]->GetLargestPossibleRegion().GetSize();
	InputImageType::IndexType imageindex;
	imageindex.Fill(0);
	InputImageType::RegionType region;
	region.SetSize(imagesize);
	region.SetIndex(imageindex);

	double max_values[MAX_CHANNELS];
	double mean_values[MAX_CHANNELS];
	int mean_values_count = 0;
	for(int counter=0; counter<n; counter++)
	{
		if(DO_FILTERING==1)
		{
				radius[0]=1;
	radius[1]=1;
	radius[2]=0;
		printf("\tPerforming median filtering on channel %d ...",counter+1);
		filt[counter]=MedianFilterType::New();
		filt[counter]->SetRadius(radius);
		filt[counter]->SetInput(im[counter]);
		filt[counter]->Update();
		iterator[counter]=IteratorType(filt[counter]->GetOutput(),filt[counter]->GetOutput()->GetLargestPossibleRegion());
		}
		else if(DO_FILTERING ==2)
		{
			printf("\tPerforming mean filtering on channel %d ...",counter+1);
			radius[0] = 2;
			radius[1] = 2;
			radius[2] = 1;
			
			filt1[counter] = MeanFilterType::New();
			filt1[counter]->SetInput(im[counter]);
			filt1[counter]->SetRadius(radius);
			filt1[counter]->Update();
			iterator[counter]=IteratorType(filt1[counter]->GetOutput(),filt1[counter]->GetOutput()->GetLargestPossibleRegion());
		}
		else if(DO_FILTERING == 0)
		{
		iterator[counter] = IteratorType(im[counter],im[counter]->GetLargestPossibleRegion());
		}

		iterator[counter].GoToBegin();
		max_values[counter]=-1;
		printf(" Done.\n",counter+1);
		for(;!iterator[counter].IsAtEnd();++iterator[counter])
		{
			if(max_values[counter]<iterator[counter].Value())
				max_values[counter] = iterator[counter].Value();
			mean_values[counter] += iterator[counter].Value();
			if(counter==0)
				mean_values_count++;
		}
		//	printf("Max%d = %lf\n",counter,max_values[counter]);
		iterator[counter].GoToBegin();
	}

	for(int counter = 0; counter < n; counter++)
	{
		char buff1[1024];
		sprintf(buff1, "C:\\Users\\arun\\Research\\unmixing_work\\build\\release\\debug_ch%d.tif",counter+1);
		//writeImage<InputImageType>(im[counter],buff1);
	}
	for(int counter = 0; counter < n; counter++)
	{
		mean_values[counter] = 0;///=mean_values_count;
		printf("mean_values[%d] = %lf\n",counter,mean_values[counter]);
	}


	// lets choose random 100000 elements to find the clusters


	vnl_matrix<double> mixed(NUM_ELEMENTS_FOR_UNMIXING,n);


	int pc = 0;
	int npc = 0;
	int total_voxels = size[0]*size[1]*size[2];

	printf("Just before matrix creation for clustering\n");
	for(;!iterator[0].IsAtEnd();)
	{

		if(rand()*1.0/RAND_MAX< 2*NUM_ELEMENTS_FOR_UNMIXING*1.0/total_voxels)
		{
			double norm = 0;
			double inf_norm = 0;
			int counter;
			for(counter=0; counter<n; counter++)
			{
				norm += iterator[counter].Get()*1.0*iterator[counter].Get();
				inf_norm = MAX(inf_norm,iterator[counter].Get());
#ifdef NO_SATURATED_PIXELS_IN_UNMIXIMG
				if(iterator[counter].Get()==255)
					break;
#endif
			}
			if(counter!=n)
			{
				for(int counter1=0; counter1< n; counter1++)
				{
					++iterator[counter1];
				}
				continue;
			}
			norm = sqrt(norm);
			//if(norm < MIN_NORM)
			if(inf_norm < MIN_NORM)
			{
				for(int counter1=0; counter1< n; counter1++)
				{
					++iterator[counter1];
				}
				continue;
			}
			for(int counter=0; counter < n; counter++)
				mixed[pc][counter] = iterator[counter].Get()-mean_values[counter];
			pc++;
		//	printf("I came into the rand\n");
		}
		else
		{
			npc++;
		//	printf("I'm in else of rand()\n");
		}
		for(int counter=0; counter< n; counter++)
		{
			++iterator[counter];
		}
		if(pc > NUM_ELEMENTS_FOR_UNMIXING-1)
			break;
	}

	printf( "npc = %d pc = %d npc+pc = %d",npc,pc,npc+pc);
	printf("Matrix created. pc = %d\n",pc);

	mixed = mixed.extract(pc,n);

	for(int counter=0; counter< n; counter++)
		iterator[counter].GoToBegin();

	vnl_matrix<double> start(n,m);

	for(int counter=0; counter<m; counter++)
	{
		for(int ncounter = 0; ncounter < n; ncounter++)
		{
			if(counter!=ncounter)
				start[ncounter][counter] = rand();
			else
				start[ncounter][counter] = rand();
		}
	}

	start.print(std::cout);

	printf("Just before unmix_cluster_from_matrix\n");
	vnl_vector <unsigned char> indices(mixed.rows()); 
	// takes a matrix and an initial estimate of eigenvectors in start. Returns the final eigenvectors in start
	//unmix_cluster_from_matrix_v2(mixed,start,indices);
#ifdef GENERATE_FINGERPRINT_MODE
	//FILE *fpp = fopen("L:\\unmixing\\SinglePositives\\Spectral\\spectral_fingerprints.txt","a+");
	FILE *fpp = fopen("L:\\unmixing\\SinglePositives\\multitrack\\spectral_fingerprints.txt","a+");
	for(int counter = 0; counter < start.columns(); counter++)
	{
		for(int counter1 = 0; counter1 < start.rows(); counter1++)
		{
			fprintf(fpp,"%lf ",start[counter1][counter]);
		}
		fprintf(fpp,"\n");
	}
	fclose(fpp);

	_exit(0);
#endif
	start = getFingerPrintMatrix();
	IteratorType oiter[20];

	for(int counter=0; counter<m; counter++)
	{
		om[counter]= InputImageType::New();
		om[counter]->SetRegions(im[0]->GetLargestPossibleRegion());
		om[counter]->Allocate();
		oiter[counter] = IteratorType(om[counter],om[counter]->GetLargestPossibleRegion());
		oiter[counter].GoToBegin();
	}
	InputImageType::Pointer assignment_image = InputImageType::New();
	assignment_image->SetRegions(om[0]->GetLargestPossibleRegion());
	assignment_image->Allocate();

	IteratorType assignment_iter = IteratorType(assignment_image,assignment_image->GetLargestPossibleRegion());
	assignment_iter.GoToBegin();

	int num_processed = 0;

	printf("Performing unmixing...\n");
	start.print(std::cout);
	//scanf("%*d");

	/*for(int counter = 0; counter < start.rows(); counter++)
	{
		double maxval = -10;
		double maxpos = -1;
		for(int counter1 = 0; counter1 < start.columns(); counter1++)
		{
			if(start[counter][counter1] > maxval)
			{
				maxpos = counter1;
				maxval = start[counter][counter1];
			}
		}
			if(maxpos != counter)
			{
				vnl_vector<double> vclvec1 = start.get_column(maxpos);
				vnl_vector<double> vclvec2 = start.get_column(counter);
				start.set_column(maxpos,vclvec2);
				start.set_column(counter,vclvec1);
			}
	}*/
printf("\n");
start.normalize_columns();
//my_normalize_columns(start);
	start.print(std::cout);
	int final_counter=0;
	vnl_matrix<double> inveigen = vnl_matrix_inverse<double>(start).pinverse(6);
	//inveigen.normalize_columns();
	inveigen.print(std::cout);
#ifdef NORMALIZE_OUTPUT_ROWS 
	printf("after normalizing rows..\n");
	inveigen.normalize_rows();
	//my_normalize_rows(inveigen);
	inveigen.print(std::cout);
#endif
	vnl_qr<double> qrdecomp(start);
	

	vnl_vector<double> inputmixed(n);
//#define UNMIX_CLUSTERING_PIXEL_CLASSIFICATION

	double min_values[MAX_CHANNELS];
	for(int counter =0; counter < m; counter++)
		min_values[counter] = 1000000.0f;

	
	for(int counter = 0; counter <n ; counter++)
	{
		if(DO_FILTERING == 0)
		{
			iterator[counter] = IteratorType(im[counter],im[counter]->GetLargestPossibleRegion());
		}
		else if (DO_FILTERING == 1)
		{
			iterator[counter] = IteratorType(filt[counter]->GetOutput(),filt[counter]->GetOutput()->GetLargestPossibleRegion());
		}
		else
		{
			iterator[counter] = IteratorType(filt1[counter]->GetOutput(),filt1[counter]->GetOutput()->GetLargestPossibleRegion());
		}
		iterator[counter].GoToBegin();
	}
	int negcount = 0;
	for(;!iterator[0].IsAtEnd();)
	{
		

		double max_proj = -123213;
		int maxpos = -1;
#if defined(UNMIX_CLUSTERING_PIXEL_CLASSIFICATION) //|| defined(UNMIX_METHOD3)
		for(int co = 0; co < m; co++)
		{
			double proj_sum = 0;
			for(int co1 = 0; co1 < n; co1++)
			{
				proj_sum += start[co1][co]*iterator[co1].Get();
			}
			if(max_proj < proj_sum)
			{
				max_proj = proj_sum;
				maxpos = co;
			}
		}
		//printf("loop 1 complete\n");
		assignment_iter.Set(maxpos);
		++assignment_iter;

#endif
		/*
		for(int co =0; co<n; co++)
		{
			inputmixed[co] = iterator[co].Get();
		}
		vnl_vector<double> outputmixed = inputmixed*inveigen;
		for(int co = 0; co < m; co++)
		{
			min_values[co] = MIN(min_values[co],outputmixed[co]);
		}*/
		for(int counter =0; counter < n; counter++)
		{
			++iterator[counter];
		}



	}

	for(int counter = 0; counter < n; counter++)
		iterator[counter].GoToBegin();
	for(int counter = 0; counter < m; counter++)
		printf("min_values[%d] = %lf\n",counter, min_values[counter]);
#ifdef UNMIX_CLUSTERING_PIXEL_CLASSIFICATION
	typedef itk::MedianImageFilter<InputImageType,InputImageType> GrayScaleMedianFilterType;
		GrayScaleMedianFilterType::Pointer binmedfilt = GrayScaleMedianFilterType::New();
		binmedfilt->SetInput(assignment_image);
		radius[0] = 1;
		radius[1] = 1;
		radius[2] = 0;
		binmedfilt->SetRadius(radius);
		binmedfilt->Update();
		assignment_image = binmedfilt->GetOutput();
		assignment_iter = IteratorType(assignment_image,assignment_image->GetLargestPossibleRegion());
		assignment_iter.GoToBegin();
		//writeImage<InputImageType>(assignment_image,"C:\\Users\\arun\\Research\\unmixing_work\\build\\release\\debug_assign_image.tif");
#endif

#ifdef UNMIX_METHOD4
/*	IloEnv env;
	IloObjective obj = IloMaximize(env);
	IloNumVarArray x(env,m,0,255);
	IloRangeArray c(env);
	for(int counter = 0; counter < n; counter++)
	{
		c.add(IloRange(env,0,1));
		for(int counter1 = 0; counter1 < m; counter1++)
		{
			c[counter].setLinearCoef(x[counter1],start[counter][counter1]);
		}
	}
	IloModel model(env);
	model.add(obj);
	model.add(c);
	IloNumArray vals(env);*/
		inveigen = start.transpose()*vnl_matrix_inverse<double>(start*start.transpose());
		inveigen.print(std::cout);
		//inveigen.normalize_rows();
#endif
	vnl_vector<int> class_count(m,0);
	for(;!oiter[0].IsAtEnd();)
	{
		/*
		   if(num_processed % 100 ==0)
		   printf("Processed %0.2lf%% voxels\r",100.0/total_voxels*num_processed);
		   */
		//		if(100.0/total_voxels*num_processed> 100)
		//			break;
		
		final_counter++;
		if(final_counter%1000000==0)
			printf("%.2lf\r",final_counter*100.0/total_voxels);

#if defined(UNMIX_CLUSTERING_PIXEL_CLASSIFICATION)
		double max_proj = -133434;
		int maxpos = -1;
		double proj_sum = 0;
		for(int co1 = 0; co1 < n; co1++)
		{
			proj_sum += start[co1][assignment_iter.Get()]*iterator[co1].Get();
		}
		;
		//printf("loop 1 complete\n");
		maxpos = assignment_iter.Get();
		//if(maxpos!=0)
		//	printf("Its not zero %d\n", rand());
		for(int co = 0; co < m; co++)
		{
			if(maxpos != co)
			{
				oiter[co].Set(0);
			}
			else
			{
				//printf("I did come here\n");
				if(proj_sum<255)
					oiter[co].Set(proj_sum);
				else
					oiter[co].Set(255);
			}
		}
		
		++assignment_iter;
#elif defined(UNMIX_METHOD3)
		
		vnl_vector<double> b(n);
		vnl_vector<int> marker(m,0);
		for(int co_n = 0; co_n < n; co_n++)
		{
			b[co_n] = iterator[co_n].Get();
		}
		for(int co_n = 0; co_n < n; co_n++)
		{
			int max_proj_co = -1;
			double max_proj_sum = 1;
			for(int co_m = 0; co_m < m ; co_m++)
			{
				double proj_sum = 0;

				for(int co1 = 0; co1 < n; co1++)
				{
					proj_sum += start[co1][co_m]*b[co1];
				}
				if(proj_sum > max_proj_sum)
				{
					max_proj_sum = proj_sum;
					max_proj_co = co_m;
				}
			}
			if(max_proj_co >=0)
			{
				//if(marker[max_proj_co]==0)
				{
					oiter[max_proj_co].Set(MAX(MIN(255,max_proj_sum),0));
					++class_count[max_proj_co];
					marker[max_proj_co] = 1;
					b = b - max_proj_sum*start.get_column(max_proj_co);
				}
				/*
				else
				{
					printf("Why am I setting it again?\n");
				}*/
			}
			else
				break;
		}
#elif defined(UNMIX_METHOD4)

		vnl_vector<double> b(n);
		for(int co =0; co<n; co++)
		{
			b[co] = iterator[co].Get()-mean_values[co];
		}
		int max_proj_co = -1;
		double max_proj_sum = 1;
		for(int co_m = 0; co_m < m ; co_m++)
		{
			double proj_sum = 0;
			for(int co1 = 0; co1 < n; co1++)
			{
				proj_sum += start[co1][co_m]*b[co1];
			}
			if(proj_sum > max_proj_sum)
			{
				max_proj_sum = proj_sum;
				max_proj_co = co_m;
			}
		}
		if(max_proj_co!=-1)
		{
			vnl_vector<double> sol(m,0);
			vnl_vector<double> solold(m,0);
			sol[max_proj_co] = b.two_norm();

			bool converged = false;
			vnl_vector<double> r(n,0);
			vnl_vector<double> d(m,0);
			int iterations;
			for(iterations = 0; iterations < 200; iterations++)
			{
				//printf("hi\n");
				r = b - start*sol;
				//printf("hi1\n");
				d = inveigen*r;
				//printf("hi2\n");
				
				sol = sol + d;
				//std::cout<<sol<<std::endl;
				//printf("hi3\n");
				for(int counter = 0; counter < sol.size(); counter++)
				{
					if(sol[counter]<0)
						sol[counter]=0;
				}
				if((sol-solold).inf_norm()<1)
				{
					//printf("hi4\n");
					break;
				}
				solold = sol;
				//printf("hi5\n");
			}
		//	if(iterations == 100)
		//		printf("hit limit %d\n",rand());
			//printf("finished\n");
			for(int co_m = 0; co_m < m; co_m++)
				oiter[co_m].Set(MAX(MIN(sol[co_m]/3,255),0));
		}
#else
		for(int co =0; co<n; co++)
		{
			inputmixed[co] = iterator[co].Get()-mean_values[co];
		}
		vnl_vector<double> outputmixed = inveigen*inputmixed;
		//vnl_vector<double> outputmixed = qrdecomp.solve(inputmixed);
		double negflag = false;
		for(int co = 0; co < m; co++)
		{	
			if(outputmixed[co]<0)
				negflag = true;
			oiter[co].Set(MAX(MIN((outputmixed[co]),255),0));
		}
		if(negflag)
			negcount++;
#endif
		
		for(int co = 0; co<n; co++)
		{
			++iterator[co];
		}
		for(int co = 0; co<m; co++)
		{
			++oiter[co];
		}
	}

	for(int count = 0; count < m; count ++)
	{
		printf("class_count[%d] = %d\n", count,class_count[count]);
	}

	/*InputImageType::Pointer processed[MAX_CHANNELS];
	for(int counter =0; counter < n; counter++)
	{
		processed[counter] = filt[counter]->GetOutput();
	}*/
	//Output2DImageType::Pointer strip1 = getGalleryView(processed,1,MAX(n,m));
	
	//writeImage<Output2DImageType>(strip1,"InputProcessedGallery.bmp");

	printf("negcount = %d\n",negcount);
	printf(" Done.\n");
	printf("End of function\n");
}
void unmix_median(InputImageType::Pointer im[],InputImageType::Pointer om[],int n)
{
	printf("Unnmixing %d channels ...\n",n);
	InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();

	printf("I'm here\n");
	MedianFilterType::Pointer filt[15];
	IteratorType iterator[15];
	IteratorType assigniter[15];

	InputImageType::SizeType radius;
	radius[0]=1;
	radius[1]=1;
	radius[2]=1;

	InputImageType::SizeType imagesize = im[0]->GetLargestPossibleRegion().GetSize();
	InputImageType::IndexType imageindex;
	imageindex.Fill(0);
	InputImageType::RegionType region;
	region.SetSize(imagesize);
	region.SetIndex(imageindex);

	double max_values[15];
	for(int counter=0; counter<n; counter++)
	{
		printf("\tPerforming median filtering on channel %d ...",counter+1);
		filt[counter]=MedianFilterType::New();
		filt[counter]->SetRadius(radius);
		filt[counter]->SetInput(im[counter]);
		filt[counter]->Update();
		om[counter]=filt[counter]->GetOutput();
		iterator[counter]=IteratorType(om[counter],om[counter]->GetLargestPossibleRegion());
		iterator[counter].GoToBegin();
		max_values[counter]=-1;
		printf(" Done.\n",counter+1);
		for(;!iterator[counter].IsAtEnd();++iterator[counter])
		{
			if(max_values[counter]<iterator[counter].Value())
				max_values[counter] = iterator[counter].Value();

		}
		//	printf("Max%d = %lf\n",counter,max_values[counter]);
		iterator[counter].GoToBegin();
	}





	int total_voxels = size[0]*size[1]*size[2];
	int num_processed = 0;
	printf("\tComputing maximum among channels ... ");
	for(;!iterator[0].IsAtEnd();)
	{
		num_processed++;
		/*
		   if(num_processed % 100 ==0)
		   printf("Processed %0.2lf%% voxels\r",100.0/total_voxels*num_processed);
		   */
		//		if(100.0/total_voxels*num_processed> 100)
		//			break;
		double max = -1;
		int maxpos = -1;
		for(int co = 0; co < n; co++)
		{
			double temp = iterator[co].Value();///max_values[co];
			if(max < temp)
			{
				max = temp;
				maxpos = co;
			}
		}
		for(int co = 0; co < n; co++)
		{
			if(maxpos != co)
			{
				iterator[co].Set(0);
			}
		}
		for(int co = 0; co<n; co++)
		{
			++iterator[co];

		}
	}
	printf(" Done.\n");
}





// main
int main(int argc, char** argv)
{
	clock_t start, end;
	start = clock();
	int status = EXIT_SUCCESS;
	if(argc< 5)
	{
		printf("Usage: unmix n m inputfile1 inputfile2 [inputfile3 ....] outputfile1 outputfile2 [outputfile3 ...]\n");
		return 1;
	}

	int n = atoi(argv[1]);
	int m = atoi(argv[2]);

	InputImageType::Pointer im[MAX_CHANNELS];
	InputImageType::Pointer om[MAX_CHANNELS];
	InputImageType::Pointer assign[MAX_CHANNELS];

	for(int counter=1; counter < n+1; counter++)
	{
		im[counter-1] = readImage<InputImageType>(argv[counter+2]);
	}

	printf("Initiating unmixing\n");
	unmix_clustering(im,om,n,m);
	printf("Finished unmix_clustering\n");
	printf("About to write output files to disk..\n");

	for(int counter=1; counter<m+1; counter++)
	{
		status |= writeImage<OutputImageType>(om[counter-1],argv[n+counter+2]);
	}

	// making the two list of images to be of equal sizes by adding empty images
	if(n<m)
	{
		for(int counter = n; counter < m; counter++)
		{
			im[counter] = getEmpty<InputImageType>(im[0]);
		}
	}
	else if (m < n)
	{
		for(int counter = m; counter < n; counter++)
		{
			om[counter] = getEmpty<InputImageType>(om[0]);
		}
	}

	Output2DImageType::Pointer strips[2];
	strips[0] = getGalleryView(im,1,MAX(n,m));
	strips[1] = getGalleryView(om,1,MAX(n,m));
	writeImage<Output2DImageType>(strips[0],"InputGallery.bmp");
	writeImage<Output2DImageType>(strips[1],"OutputGallery.bmp");
	return 0;

}
