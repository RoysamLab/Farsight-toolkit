#include "ftkSpectralUnmixing.h"

namespace ftk
{
double square_function (double a)
{
	return a*a;
}

//Constructor
SpectralUnmixing::SpectralUnmixing()
{
	UnmixedImage = ftk::Image::New();
	this->unmixMode = PIXEL_CLASSIFICATION;
}
SpectralUnmixing::~SpectralUnmixing()
{
}
void SpectralUnmixing::SetInputImage(ftk::Image::Pointer image)
{
	this->Image = image;
	this->NChannels = (int)Image->GetImageInfo()->numChannels;
}
void SpectralUnmixing::SetNumberOfChannels(int m)
{
	this->MChannels = m;
}
void SpectralUnmixing::SetUnmixMode(UnmixMode mode)
{
	this->unmixMode = mode;
	if(unmixMode == LINEAR_UNMIX)
	{
		if(NChannels>=MChannels)
			this->sysMode = OVER_DETERMINED;
		else
			this->sysMode = UNDER_DETERMINED;
	}
}


void SpectralUnmixing::Update(void)
{
	if(MChannels>MAX_CHANNS||NChannels>MAX_CHANNS)
	{
		std::cout<<"Cannot handle more than: "<<MAX_CHANNS<<" channels.\n";
		return;
	}
	// Don't forget to pixel data type
	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(2); //DEEP_COPY mode

	// First Get the FingerPrint Matrix:
	//InputImageType::Pointer imdata[MAX_CHANNS];
	//for(int ch = 0; ch< NChannels; ++ch)
	//{
	//	imdata[ch] = Image->GetItkPtr<InputPixelType>(0,ch,mode);
	//}
//	vnl_matrix<double> FingerPrintMatrix = GetFingerPrintMatrix(imdata);
	vnl_matrix<double> FingerPrintMatrix = this->GetFingerPrintMatrix();
	printf("Finished Estimating from first Image, I am gonna use it to unmix everything else\n");
	
	int numTSlices = (int)Image->GetImageInfo()->numTSlices;
	for(int T=0; T<numTSlices; ++T)
	{
		InputImageType::Pointer im[MAX_CHANNS];
		InputImageType::Pointer om[MAX_CHANNS];

		for(int ch = 0; ch< NChannels; ++ch)
		{
			im[ch] = Image->GetItkPtr<InputPixelType>(T,ch,mode);
		}
		printf("Initiating unmixing for T = %d\n",T);
		UnmixClustering(im,om,FingerPrintMatrix);
		printf("Finished UnmixClustering for T = %d\n",T);
	}
	this->ConvertOutputToftk();
	

}
//vnl_matrix<double> SpectralUnmixing::GetFingerPrintMatrix(InputImageType::Pointer im[])
vnl_matrix<double> SpectralUnmixing::GetFingerPrintMatrix(void)
{
	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(2); //DEEP_COPY mode
	InputImageType::Pointer im[MAX_CHANNS];
	int numTSlices = (int)Image->GetImageInfo()->numTSlices;
	vnl_matrix<double> mixed(NUM_ELEMENTS_FOR_UNMIXING,NChannels);
	int pc = 0;
	int npc = 0;
	for(int T=0; T<numTSlices; ++T)
	{

		printf("Selecting from T:%d\n",T);
		for(int ch = 0; ch< NChannels; ++ch)
		{
			im[ch] = Image->GetItkPtr<InputPixelType>(0,ch,mode);
		}
		
		InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
		IteratorType iterator[MAX_CHANNS];
		IteratorType assigniter[MAX_CHANNS];

		double max_values[MAX_CHANNS];
		double mean_values[MAX_CHANNS];
		int mean_values_count = 0;
		// Compute max and mean values of the images:
		for(int counter=0; counter<NChannels; counter++) 
		{
			iterator[counter] = IteratorType(im[counter],im[counter]->GetLargestPossibleRegion());
			iterator[counter].GoToBegin();
			max_values[counter]=-1;
		//	printf(" Done.\n",counter+1);
			for(;!iterator[counter].IsAtEnd();++iterator[counter]) // compute max and mean value
			{
				if(max_values[counter]<iterator[counter].Value())
					max_values[counter] = iterator[counter].Value();
				mean_values[counter] += iterator[counter].Value();
				if(counter==0)
					mean_values_count++;
			}
		//	printf("Max%d = %lf\n",counter,max_values[counter]);
		//	printf("mean_values[%d] = %lf\n",counter,mean_values[counter]);
			iterator[counter].GoToBegin();
		}

	// Choose the voxels data for the unmixing algorithm( thesis: Yi)and put them in mixed(NUM_ELEMENTS_FOR_UNMIXING,n).

		int total_voxels = size[0]*size[1]*size[2];
		printf("Just before matrix creation for clustering\n");
		for(;!iterator[0].IsAtEnd();)		// iterate through all voxels
		{
			if(rand()*1.0/RAND_MAX< 2*NUM_ELEMENTS_FOR_UNMIXING*1.0/total_voxels)
			{
				double norm = 0;
				double inf_norm = 0;
				int counter;
				for(counter=0; counter<NChannels; counter++)
				{
					norm += iterator[counter].Get()*1.0*iterator[counter].Get();
					inf_norm = MAX(inf_norm,iterator[counter].Get());
				}
				if(counter!=NChannels)
				{
					for(int counter1=0; counter1< NChannels; counter1++)	// increment the iterator to the next voxel
					{
						++iterator[counter1];
					}
					continue;
				}
				norm = sqrt(norm);
				if(inf_norm < MIN_NORM)
				{
					for(int counter1=0; counter1< NChannels; counter1++)	// increment the iterator the next voxel
					{
						++iterator[counter1];
					}
					continue;
				}
				for(int counter=0; counter < NChannels; counter++)
					mixed[pc][counter] = iterator[counter].Get()-mean_values[counter];
				pc++;
				//printf("I came into the rand\n");
			}
			else
			{
				npc++;
				//printf("I'm in else of rand()\n");
			}
			for(int counter=0; counter< NChannels; counter++)
			{
				++iterator[counter];
			}
			if(pc > NUM_ELEMENTS_FOR_UNMIXING-1)
				break;
		}
		if(pc > NUM_ELEMENTS_FOR_UNMIXING-1)
			break;
	}
	// Now Compute the fingerprint matrix:
	printf("Finished Collecting Mixed Data\n");
	mixed = mixed.extract(pc,NChannels);
	printf( "npc = %d pc = %d npc+pc = %d",npc,pc,npc+pc);
	printf("Matrix created. pc = %d\n",pc);
	// Generate a matrix of random numbers:
	vnl_matrix<double> start(NChannels,MChannels);
	for(int counter=0; counter<MChannels; counter++)
	{
		for(int ncounter = 0; ncounter < NChannels; ncounter++)
		{
			start[ncounter][counter] = rand();
		}
	}
	start.print(std::cout);
	printf("Just before EstimateFingerPrintMatrix\n");
	vnl_vector <unsigned char> indices(mixed.rows()); 


	this->EstimateFingerPrintMatrix(mixed,start,indices);
	printf("Finished estimating Finger Print Matrix:\n");
	start.print(std::cout);
	start.normalize_columns();
	printf("Normalized Finger Print Matrix:\n");
	start.print(std::cout);
	return start;
}
void SpectralUnmixing::UnmixClustering(InputImageType::Pointer im[],InputImageType::Pointer om[],vnl_matrix<double> start)
{
	if(unmixMode == PIXEL_CLASSIFICATION)
		this->UnmixPureChannels(im,om,start);
	else if(unmixMode == LINEAR_UNMIX && sysMode == OVER_DETERMINED) 
		this->UnmixUsingPseudoInverse(im,om,start);
	else 
		this->UnmixUsingIterations(im,om,start);

}
 
void SpectralUnmixing::UnmixUsingPseudoInverse(InputImageType::Pointer im[],InputImageType::Pointer om[],vnl_matrix<double> &start)
{
	InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
	int total_voxels = size[0]*size[1]*size[2];
	IteratorType iterator[MAX_CHANNS];
	IteratorType assigniter[MAX_CHANNS];
	IteratorType oiter[20];
	for(int counter=0; counter<MChannels; counter++)	
	{
		om[counter]= InputImageType::New();
		om[counter]->SetRegions(im[0]->GetLargestPossibleRegion());
		om[counter]->Allocate();
		oiter[counter] = IteratorType(om[counter],om[counter]->GetLargestPossibleRegion());
		oiter[counter].GoToBegin();
	}
	// allocate memeory for the assignment image
	InputImageType::Pointer assignment_image = InputImageType::New();
	assignment_image->SetRegions(om[0]->GetLargestPossibleRegion());
	assignment_image->Allocate();
	IteratorType assignment_iter = IteratorType(assignment_image,assignment_image->GetLargestPossibleRegion());
	assignment_iter.GoToBegin();
	int num_processed = 0;
	printf("Performing unmixing...\n");
	int final_counter=0;
	vnl_matrix<double> inveigen = vnl_matrix_inverse<double>(start).pinverse(6); //get the pseudo inverse of start (a random matrix)
	inveigen.print(std::cout);
	vnl_qr<double> qrdecomp(start);	// get the QR decomposition of start( supposedly the fingerprint matrix)
	vnl_vector<double> inputmixed(NChannels);
	double min_values[MAX_CHANNS];
	for(int counter =0; counter <MChannels; counter++)
		min_values[counter] = 1000000.0f;
	for(int counter = 0; counter <NChannels ; counter++)
	{
		iterator[counter] = IteratorType(im[counter],im[counter]->GetLargestPossibleRegion());
		iterator[counter].GoToBegin();
	}
	int negcount = 0;
	for(;!iterator[0].IsAtEnd();) // iterate through voxels
	{
		double max_proj = -123213;
		int maxpos = -1;
		for(int counter =0; counter < NChannels; counter++)
		{
			++iterator[counter];
		}
	}
	for(int counter = 0; counter < NChannels; counter++)
		iterator[counter].GoToBegin();
	for(int counter = 0; counter < MChannels; counter++)
		printf("min_values[%d] = %lf\n",counter, min_values[counter]);
	vnl_vector<int> class_count(MChannels,0);
	for(;!oiter[0].IsAtEnd();)
	{
		final_counter++;
		if(final_counter%1000000==0)
			printf("%.2lf\r",final_counter*100.0/total_voxels);
		for(int co =0; co<NChannels; co++)
		{
		//	inputmixed[co] = iterator[co].Get()-mean_values[co];
			inputmixed[co] = iterator[co].Get();
		}
		vnl_vector<double> outputmixed = inveigen*inputmixed;
		//vnl_vector<double> outputmixed = qrdecomp.solve(inputmixed);
		double negflag = false;
		for(int co = 0; co < MChannels; co++)
		{	
			if(outputmixed[co]<0)
				negflag = true;
			oiter[co].Set(MAX(MIN((outputmixed[co]),255),0));
		}
		if(negflag)
			negcount++;
		for(int co = 0; co<NChannels; co++)
		{
			++iterator[co];
		}
		for(int co = 0; co<MChannels; co++)
		{
			++oiter[co];
		}
	}
	// store the pointers to the channels in a vector
	std::vector<InputImageType::Pointer> tmp;
	for(int co = 0; co<MChannels; co++)
		tmp.push_back(om[co]);
	Unmixed_Images.push_back(tmp);
	for(int count = 0; count < MChannels; count ++)
	{
		printf("class_count[%d] = %d\n", count,class_count[count]);
	}
	printf(" Done.\n");
	printf("End of function\n");
}
void SpectralUnmixing::UnmixUsingIterations(InputImageType::Pointer im[],InputImageType::Pointer om[],vnl_matrix<double> &start)
{
	InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
	int total_voxels = size[0]*size[1]*size[2];
	IteratorType iterator[MAX_CHANNS];
	IteratorType assigniter[MAX_CHANNS];
	IteratorType oiter[20];
	for(int counter=0; counter<MChannels; counter++)	
	{
		om[counter]= InputImageType::New();
		om[counter]->SetRegions(im[0]->GetLargestPossibleRegion());
		om[counter]->Allocate();
		oiter[counter] = IteratorType(om[counter],om[counter]->GetLargestPossibleRegion());
		oiter[counter].GoToBegin();
	}
	// allocate memeory for the assignment image
	InputImageType::Pointer assignment_image = InputImageType::New();
	assignment_image->SetRegions(om[0]->GetLargestPossibleRegion());
	assignment_image->Allocate();
	IteratorType assignment_iter = IteratorType(assignment_image,assignment_image->GetLargestPossibleRegion());
	assignment_iter.GoToBegin();
	int num_processed = 0;
	printf("Performing unmixing...\n");
	int final_counter=0;
	vnl_vector<double> inputmixed(NChannels);
	double min_values[MAX_CHANNS];
	for(int counter =0; counter <MChannels; counter++)
		min_values[counter] = 1000000.0f;
	for(int counter = 0; counter <NChannels ; counter++)
	{
		iterator[counter] = IteratorType(im[counter],im[counter]->GetLargestPossibleRegion());
		iterator[counter].GoToBegin();
	}

	// Unmix Using iterative method(Under-determined System):
	int negcount = 0;
	for(;!iterator[0].IsAtEnd();) // iterate through voxels
	{
		double max_proj = -123213;
		int maxpos = -1;
		for(int counter =0; counter < NChannels; counter++)
		{
			++iterator[counter];
		}
	}
	for(int counter = 0; counter < NChannels; counter++)
		iterator[counter].GoToBegin();
	for(int counter = 0; counter < MChannels; counter++)
		printf("min_values[%d] = %lf\n",counter, min_values[counter]);
	vnl_vector<int> class_count(MChannels,0);
	for(;!oiter[0].IsAtEnd();)
	{
		final_counter++;
		if(final_counter%1000000==0)
			printf("%.2lf\r",final_counter*100.0/total_voxels);
				
		vnl_vector<double> b(NChannels);
		vnl_vector<int> marker(MChannels,0);
		for(int co_n = 0; co_n < NChannels; co_n++)
		{
			b[co_n] = iterator[co_n].Get();
		}
		for(int co_n = 0; co_n < NChannels; co_n++)
		{
			int max_proj_co = -1;
			double max_proj_sum = 1;
			for(int co_m = 0; co_m < MChannels ; co_m++)
			{
				double proj_sum = 0;

				for(int co1 = 0; co1 < NChannels; co1++)
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
		for(int co = 0; co<NChannels; co++)
		{
			++iterator[co];
		}
		for(int co = 0; co<MChannels; co++)
		{
			++oiter[co];
		}
	}
	// store the pointers to the channels in a vector
	std::vector<InputImageType::Pointer> tmp;
	for(int co = 0; co<MChannels; co++)
		tmp.push_back(om[co]);
	Unmixed_Images.push_back(tmp);
	for(int count = 0; count < MChannels; count ++)
	{
		printf("class_count[%d] = %d\n", count,class_count[count]);
	}
	printf(" Done.\n");
	printf("End of function\n");
}



void SpectralUnmixing::UnmixPureChannels(InputImageType::Pointer im[],InputImageType::Pointer om[],vnl_matrix<double> &start)
{
	InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
	int total_voxels = size[0]*size[1]*size[2];
	IteratorType iterator[MAX_CHANNS];
	IteratorType assigniter[MAX_CHANNS];

	IteratorType oiter[20];
	for(int counter=0; counter<MChannels; counter++)	
	{
		om[counter]= InputImageType::New();
		om[counter]->SetRegions(im[0]->GetLargestPossibleRegion());
		om[counter]->Allocate();
		oiter[counter] = IteratorType(om[counter],om[counter]->GetLargestPossibleRegion());
		oiter[counter].GoToBegin();
	}
	// allocate memeory for the assignment image
	InputImageType::Pointer assignment_image = InputImageType::New();
	assignment_image->SetRegions(om[0]->GetLargestPossibleRegion());
	assignment_image->Allocate();
	IteratorType assignment_iter = IteratorType(assignment_image,assignment_image->GetLargestPossibleRegion());
	assignment_iter.GoToBegin();
	int num_processed = 0;
	printf("Performing unmixing...\n");
	int final_counter=0;
	vnl_vector<double> inputmixed(NChannels);
	double min_values[MAX_CHANNS];
	for(int counter =0; counter <MChannels; counter++)
		min_values[counter] = 1000000.0f;
	for(int counter = 0; counter <NChannels ; counter++)
	{
		iterator[counter] = IteratorType(im[counter],im[counter]->GetLargestPossibleRegion());
		iterator[counter].GoToBegin();
	}
	// Pixel Classification Assignment
	int negcount = 0;
	for(;!iterator[0].IsAtEnd();) // iterate through voxels
	{
		double max_proj = -123213;
		int maxpos = -1;
		for(int co = 0; co < MChannels; co++)	// iterate through the output channels
		{
			double proj_sum = 0;
			for(int co1 = 0; co1 < NChannels; co1++)	// iterate through the input channels
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
		for(int counter =0; counter < NChannels; counter++)
		{
			++iterator[counter];
		}
	}
	for(int counter = 0; counter < NChannels; counter++)
		iterator[counter].GoToBegin();
	for(int counter = 0; counter < MChannels; counter++)
		printf("min_values[%d] = %lf\n",counter, min_values[counter]);

	typedef itk::MedianImageFilter<InputImageType,InputImageType> GrayScaleMedianFilterType;
	GrayScaleMedianFilterType::Pointer binmedfilt = GrayScaleMedianFilterType::New();
	binmedfilt->SetInput(assignment_image);
	InputImageType::SizeType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;
	binmedfilt->SetRadius(radius);
	binmedfilt->Update();
	assignment_image = binmedfilt->GetOutput();
	assignment_iter = IteratorType(assignment_image,assignment_image->GetLargestPossibleRegion());
	assignment_iter.GoToBegin();

	vnl_vector<int> class_count(MChannels,0);
	for(;!oiter[0].IsAtEnd();)
	{
		final_counter++;
		if(final_counter%1000000==0)
			printf("%.2lf\r",final_counter*100.0/total_voxels);
		double max_proj = -133434;
		int maxpos = -1;
		double proj_sum = 0;
		for(int co1 = 0; co1 < NChannels; co1++)
		{
			proj_sum += start[co1][assignment_iter.Get()]*iterator[co1].Get();
		}
		maxpos = assignment_iter.Get();
		for(int co = 0; co < MChannels; co++)
		{
			if(maxpos != co)
			{
				oiter[co].Set(0);
			}
			else
			{
				if(proj_sum<255)
					oiter[co].Set(proj_sum);
				else
					oiter[co].Set(255);
			}
		}
		
		++assignment_iter;
		for(int co = 0; co<NChannels; co++)
		{
			++iterator[co];
		}
		for(int co = 0; co<MChannels; co++)
		{
			++oiter[co];
		}
	}
	// store the pointers to the channels in a vector
	std::vector<InputImageType::Pointer> tmp;
	for(int co = 0; co<MChannels; co++)
		tmp.push_back(om[co]);
	Unmixed_Images.push_back(tmp);
	for(int count = 0; count < MChannels; count ++)
	{
		printf("class_count[%d] = %d\n", count,class_count[count]);
	}
	printf(" Done.\n");
	printf("End of function\n");

}
void SpectralUnmixing::ConvertOutputToftk(void)	
{
	//ftk::Image::Pointer tmp_image = ftk::Image::New();
	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(2); //DEEP_COPY mode
	std::vector< std::vector <unsigned char> > channelColors = Image->GetImageInfo()->channelColors;
	ftk::Image::DataType dataType = Image->GetImageInfo()->dataType;
	unsigned char databpPix = Image->GetImageInfo()->bytesPerPix;
	unsigned char labelbpPix = Image->GetImageInfo()->bytesPerPix;
	unsigned short cs = Image->GetImageInfo()->numColumns;
	unsigned short rs = Image->GetImageInfo()->numRows;
	unsigned short zs = Image->GetImageInfo()->numZSlices;
	std::vector< std::vector <std::string> > FileNames = Image->GetTimeChannelFilenames();

	printf("converting to ftk (time): 0\n");
	for(int ch =0;ch<MChannels;++ch)
	{
		
		printf("converting to ftk (channel): %d\n",ch+1);
		std::stringstream ss;
		ss<<ch;
		std::string name ="channel"+ss.str();
		std::vector <unsigned char> colors;
		this->getColor(ch,&colors);
		UnmixedImage->AppendChannelFromData3D(Unmixed_Images[0][ch]->GetBufferPointer(), dataType, databpPix, cs, rs, zs, name,colors, true);
	}


	for(int t = 1;t<(int)Unmixed_Images.size();++t)
	{
		printf("converting to ftk (time): %d\n",t);
		ftk::Image::Pointer tmp_image = ftk::Image::New();

		for(int ch =0;ch<MChannels;++ch)
		{
			printf("converting to ftk (channel): %d\n",ch+1);
			std::stringstream ss;
			ss<<ch;
			std::string name = "channel"+ss.str();
			tmp_image->AppendChannelFromData3D(Unmixed_Images[t][ch]->GetBufferPointer(), dataType, databpPix, cs, rs, zs, name, channelColors.at(ch), true);
		}
		UnmixedImage->AppendImage(tmp_image,mode,true);
	}
	std::vector< std::vector <std::string> > tmp_filenames;
	for(int i = 0; i< UnmixedImage->GetImageInfo()->numTSlices; ++i)
	{
		std::vector <std::string> tmp_file;
		for(int ch = 0; ch< UnmixedImage->GetImageInfo()->numChannels; ++ch)
		{
			std::stringstream ss;
			ss<<ch;
			std::string name = ftk::GetFilePath(FileNames.at(i).at(ch))+"\\unmixed_ch"+ss.str()+"_"+ftk::GetFilenameFromFullPath(FileNames.at(i).at(ch));
			tmp_file.push_back(name);
			std::cout<<name<<std::endl;
		}
		tmp_filenames.push_back(tmp_file);
	}
	UnmixedImage->SetTimeChannelFilenames(tmp_filenames);


}



void SpectralUnmixing::EstimateFingerPrintMatrix(vnl_matrix<double> mixed, vnl_matrix<double> &start, vnl_vector<unsigned char> &indices)
{
	printf("Entered EstimateFingerPrintMatrix\n");
	int m = start.cols();
	int n = start.rows();

//	std::ofstream f;
//	f.open("C:/Users/arun/Research/unmixing_work/build/release/matrix.txt",std::ios_base::out);
	bool converge = false;

	start.normalize_columns();
	// some neat way of generating sum(k.^2,2)*ones(1,n);
	vnl_matrix<double> sqmixed(mixed);
	sqmixed.apply(square_function);		// square the mixed matrix
	vnl_matrix<double> ones(n,1);
	ones.fill(1);
	sqmixed = sqmixed*ones;
	ones = ones.transpose();
	sqmixed = sqmixed*ones;
	// end of neat way
	vnl_vector<double> max_along_rows(mixed.rows()); // define a vector with a size of the number of samples chosen
	double minumdoublevalue = -std::numeric_limits<double>::max();

	std::vector< vnl_vector<double> > subsets[10];	// define an array of 10 std::vectors of vnl_vectors (index of the standard vector is the index of the cluster
	int iteration_number = 0;
	vnl_matrix<double> old_start;					
	vnl_matrix<double> A[10];						
	start.print(std::cout);
	//printf("Entered the loop\n");
	while(!converge)
	{
	//	printf("Finding maximum along components\n");
		iteration_number++;
		if(iteration_number>50)
			break;
		vnl_matrix<double> out = mixed*start;

		max_along_rows.fill(minumdoublevalue);
		int num_rows = mixed.rows();
		int num_columns = m;
		old_start = start;
		start.fill(0);// set the current estimate to a zero matrix.
		for(int counter=0; counter< m; counter++)
		{
			subsets[counter].clear();
		}

		for(int counterx = 0; counterx < num_rows; counterx++) // iterate through sample values of mixed
		{
			
			indices[counterx]=-1;								// indices is a vector
			for(int countery = 0; countery < num_columns; countery++) // iterate through the number of output channels
			{


				if(max_along_rows[counterx] < out[counterx][countery])
				{
					max_along_rows[counterx] = out[counterx][countery]; // max_along_rows will contain the maximum value of the projection (dot product) among the cluster vectors
					indices[counterx] = countery;						// indices will contain the index of the cluster the current vector belongs to
				}
			}
			double temp=0;
	
			subsets[indices[counterx]].push_back(mixed.get_row(counterx)); //pushback the vectors for each class 

		}
		
		for(int counter = 0; counter < m; counter++) // iterate through the number of clusters
		{
			vnl_vector<double> p(n);
			p.fill(1);
			vnl_vector<double> t(n),old_p(n);
			t.fill(0);
			// create vectors of channel size
			old_p = p;
			while(1)
			{
				t.fill(0);
				for(int counter1 = 0; counter1 < subsets[counter].size(); counter1++) // go to each cluster and iterate through the vectors in the cluster
				{
					t = t + (dot_product(subsets[counter][counter1],p)*subsets[counter][counter1]);		// sum the elements of the vector, multiply by the vector and add to t
				}
				p = t.normalize();
				if((p-old_p).two_norm()< 1e-6)
					break;
				old_p = p;
			}
			start.set_column(counter,p); // computing the cluster centers
		}
		start.normalize_columns();
		//start.print(std::cout);
		if(start.get_column(0).is_zero() == 1 || start.get_column(1).is_zero() == 1)
		{
			srand(time(NULL));
			//printf("resetting...\n");
			for(int counter=0; counter<m; counter++)
			{
				for(int ncounter = 0; ncounter < n; ncounter++)
				{
					start[ncounter][counter] = rand()%100;
				}
			}
			std::cout<<std::setprecision(2);
			//start.print(std::cout);
			start.normalize_columns();
			continue;
		}

	 	double change = (old_start-start).apply(square_function).absolute_value_sum();
		if(change < 0.00001)
			converge = true;
	}
}

void SpectralUnmixing::getColor(int numChann,std::vector<unsigned char> *channelColors)
{

	switch(numChann)
	{
		//Cyan
	case 0:
		(*channelColors).push_back(0);
		(*channelColors).push_back(255);
		(*channelColors).push_back(255);

		break;
		//Green
	case 1:
		(*channelColors).push_back(0);
		(*channelColors).push_back(255);
		(*channelColors).push_back(0);

		break;
		//Red
	case 2:
		(*channelColors).push_back(255);
		(*channelColors).push_back(0);
		(*channelColors).push_back(0);

		break;
		//Blue
	case 3:
		(*channelColors).push_back(0);
		(*channelColors).push_back(0);
		(*channelColors).push_back(255);

		break;

		//Orange 	255-165-0
	case 4:
		(*channelColors).push_back(255);
		(*channelColors).push_back(165);
		(*channelColors).push_back(0);

		break;
		//Violet 	238-130-238
	case 5:
		(*channelColors).push_back(238);
		(*channelColors).push_back(130);
		(*channelColors).push_back(238);

		break;
		//Yellow 	255-255-0
	case 6:
		(*channelColors).push_back(255);
		(*channelColors).push_back(255);
		(*channelColors).push_back(0);

		break;
		//Dark Green 	0-100-0
	case 7:
		(*channelColors).push_back(0);
		(*channelColors).push_back(100);
		(*channelColors).push_back(0);

		break;
		//Royal Blue 	65-105-225
	case 8:
		(*channelColors).push_back(65);
		(*channelColors).push_back(105);
		(*channelColors).push_back(225);

		break;
		//Gray 	190-190-190
	case 9:
		(*channelColors).push_back(190);
		(*channelColors).push_back(190);
		(*channelColors).push_back(190);

		break;
	}
}
}// end of namespace