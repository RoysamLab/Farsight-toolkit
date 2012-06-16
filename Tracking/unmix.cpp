#include "helpers.h"

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

int num_files;



struct FileData{
	std::string id;
	std::string filename;
	int channel;
	int time_point;
	void Print()
	{
		//printf("id: |%s| filename: |%s| channel: %d time_point: %d\n",id.c_str(),filename.c_str(),channel,time_point);
	}
};

#define FILES_MAX 6000
FileData files[FILES_MAX];

void unmix_median(InputImageType::Pointer im[],InputImageType::Pointer om[],int n)
{
	printf("Unnmixing %d channels ...\n",n);
	InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();

	printf("I'm here\n");
	MedianFilterType::Pointer filt[15];
	IteratorType iterator[15];
	//IteratorType assigniter[15];

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
		printf(" Done %d.\n",counter+1);
		for(;!iterator[counter].IsAtEnd();++iterator[counter])
		{
			if(max_values[counter]<iterator[counter].Value())
				max_values[counter] = iterator[counter].Value();

		}
		//	printf("Max%d = %lf\n",counter,max_values[counter]);
		iterator[counter].GoToBegin();
	}





	//int total_voxels = size[0]*size[1]*size[2];
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

void unmixMPIInternal(InputImageType::Pointer im [4],InputImageType::Pointer om[4], InputImageType::Pointer assignment[4])
{
	printf("Performing unmixing ...\n");
	InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();


	IteratorType iterator[4];
	IteratorType assigniter[4];

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


	for(int counter=0; counter<4; counter++)
	{
		printf("\tPerforming median filtering on channel %d ...",counter+1);
		om[counter]=im[counter];
		assignment[counter]=InputImageType::New();
		assignment[counter]->SetRegions(region);
		assignment[counter]->Allocate();
		assigniter[counter]=IteratorType(assignment[counter],assignment[counter]->GetLargestPossibleRegion());
		assigniter[counter].GoToBegin();
		iterator[counter]=IteratorType(om[counter],om[counter]->GetLargestPossibleRegion());
		iterator[counter].GoToBegin();
		printf(" Done %d.\n",counter+1);
	}

	//int total_voxels = size[0]*size[1]*size[2];
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
		for(int co = 0; co < 4; co++)
		{
			unsigned char temp = iterator[co].Value();
			if(max < temp)
			{
				max = temp;
				maxpos = co;
			}
		}
		for(int co = 0; co < 4; co++)
		{
			if(maxpos != co)
			{
				iterator[co].Set(0);
				assigniter[co].Set(0);
			}
			else
			{
				assigniter[co].Set(255);
			}
		}
		for(int co = 0; co<4; co++)
		{
			++iterator[co];
			++assigniter[co];
		}
	}
	printf(" Done.\n");
}

#define mxIsFinite(a) ((a)<1e6)

void assignmentsuboptimal1(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns)
{
	bool infiniteValueFound, finiteValueFound, repeatSteps, allSinglyValidated, singleValidationFound;
	int n, row, col, tmpRow, tmpCol, nOfElements;
	int *nOfValidObservations, *nOfValidTracks;
	double value, minValue, *distMatrix, inf;

	inf = 1e10;

	/* make working copy of distance Matrix */
	nOfElements   = nOfRows * nOfColumns;
	distMatrix    = (double *)malloc(nOfElements * sizeof(double));
	for(n=0; n<nOfElements; n++)
		distMatrix[n] = distMatrixIn[n];

	/* initialization */
	*cost = 0;
#ifdef ONE_INDEXING
	for(row=0; row<nOfRows; row++)
		assignment[row] =  0.0;
#else
	for(row=0; row<nOfRows; row++)
		assignment[row] = -1.0;
#endif

	/* allocate memory */
	nOfValidObservations  = (int *)calloc(nOfRows,    sizeof(int));
	nOfValidTracks        = (int *)calloc(nOfColumns, sizeof(int));

	/* compute number of validations */
	infiniteValueFound = false;
	finiteValueFound  = false;
	for(row=0; row<nOfRows; row++)
		for(col=0; col<nOfColumns; col++)
			if(mxIsFinite(distMatrix[row + nOfRows*col]))
			{
				nOfValidTracks[col]       += 1;
				nOfValidObservations[row] += 1;
				finiteValueFound = true;
			}
			else
				infiniteValueFound = true;

	if(infiniteValueFound)
	{
		if(!finiteValueFound)
			return;

		repeatSteps = true;

		while(repeatSteps)
		{
			repeatSteps = false;

			/* step 1: reject assignments of multiply validated tracks to singly validated observations		 */
			for(col=0; col<nOfColumns; col++)
			{
				singleValidationFound = false;
				for(row=0; row<nOfRows; row++)
					if(mxIsFinite(distMatrix[row + nOfRows*col]) && (nOfValidObservations[row] == 1))
					{
						singleValidationFound = true;
						break;
					}

					if(singleValidationFound)
					{
						for(row=0; row<nOfRows; row++)
							if((nOfValidObservations[row] > 1) && mxIsFinite(distMatrix[row + nOfRows*col]))
							{
								distMatrix[row + nOfRows*col] = inf;
								nOfValidObservations[row] -= 1;							
								nOfValidTracks[col]       -= 1;	
								repeatSteps = true;				
							}
					}
			}

			/* step 2: reject assignments of multiply validated observations to singly validated tracks */
			if(nOfColumns > 1)			
			{	
				for(row=0; row<nOfRows; row++)
				{
					singleValidationFound = false;
					for(col=0; col<nOfColumns; col++)
						if(mxIsFinite(distMatrix[row + nOfRows*col]) && (nOfValidTracks[col] == 1))
						{
							singleValidationFound = true;
							break;
						}

						if(singleValidationFound)
						{
							for(col=0; col<nOfColumns; col++)
								if((nOfValidTracks[col] > 1) && mxIsFinite(distMatrix[row + nOfRows*col]))
								{
									distMatrix[row + nOfRows*col] = inf;
									nOfValidObservations[row] -= 1;
									nOfValidTracks[col]       -= 1;
									repeatSteps = true;								
								}
						}
				}
			}
		} /* while(repeatSteps) */

		/* for each multiply validated track that validates only with singly validated  */
		/* observations, choose the observation with minimum distance */
		for(row=0; row<nOfRows; row++)
		{
			if(nOfValidObservations[row] > 1)
			{
				allSinglyValidated = true;
				minValue = inf;
				for(col=0; col<nOfColumns; col++)
				{
					value = distMatrix[row + nOfRows*col];
					if(mxIsFinite(value))
					{
						if(nOfValidTracks[col] > 1)
						{
							allSinglyValidated = false;
							break;
						}
						else if((nOfValidTracks[col] == 1) && (value < minValue))
						{
							tmpCol   = col;
							minValue = value;
						}
					}
				}

				if(allSinglyValidated)
				{
#ifdef ONE_INDEXING
					assignment[row] = tmpCol + 1;
#else
					assignment[row] = tmpCol;
#endif
					*cost += minValue;
					for(n=0; n<nOfRows; n++)
						distMatrix[n + nOfRows*tmpCol] = inf;
					for(n=0; n<nOfColumns; n++)
						distMatrix[row + nOfRows*n] = inf;
				}
			}
		}

		/* for each multiply validated observation that validates only with singly validated  */
		/* track, choose the track with minimum distance */
		for(col=0; col<nOfColumns; col++)
		{
			if(nOfValidTracks[col] > 1)
			{
				allSinglyValidated = true;
				minValue = inf;
				for(row=0; row<nOfRows; row++)
				{
					value = distMatrix[row + nOfRows*col];
					if(mxIsFinite(value))
					{
						if(nOfValidObservations[row] > 1)
						{
							allSinglyValidated = false;
							break;
						}
						else if((nOfValidObservations[row] == 1) && (value < minValue))
						{
							tmpRow   = row;
							minValue = value;
						}
					}
				}

				if(allSinglyValidated)
				{
#ifdef ONE_INDEXING
					assignment[tmpRow] = col + 1;
#else
					assignment[tmpRow] = col;
#endif
					*cost += minValue;
					for(n=0; n<nOfRows; n++)
						distMatrix[n + nOfRows*col] = inf;
					for(n=0; n<nOfColumns; n++)
						distMatrix[tmpRow + nOfRows*n] = inf;
				}
			}
		}	
	} /* if(infiniteValueFound) */


	/* now, recursively search for the minimum element and do the assignment */
	while(true)
	{
		/* find minimum distance observation-to-track pair */
		minValue = inf;
		for(row=0; row<nOfRows; row++)
			for(col=0; col<nOfColumns; col++)
			{
				value = distMatrix[row + nOfRows*col];
				if(mxIsFinite(value) && (value < minValue))
				{
					minValue = value;
					tmpRow   = row;
					tmpCol   = col;
				}
			}

			if(mxIsFinite(minValue))
			{
#ifdef ONE_INDEXING
				assignment[tmpRow] = tmpCol+ 1;
#else
				assignment[tmpRow] = tmpCol;
#endif
				*cost += minValue;
				for(n=0; n<nOfRows; n++)
					distMatrix[n + nOfRows*tmpCol] = inf;
				for(n=0; n<nOfColumns; n++)
					distMatrix[tmpRow + nOfRows*n] = inf;			
			}
			else
				break;

	} /* while(true) */

	/* free allocated memory */
	free(nOfValidObservations);
	free(nOfValidTracks);


}

//void assignmentsuboptimal1(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns)
#define USE_VNL_HUNGARIAN 
//vcl_vector<unsigned int> getTimeAssociations(std::vector<FeaturesType> &a,std::vector<FeaturesType> &b)
//{
//
//#define DEBUG_RANK -1
//	int rows =0;
//	int cols =0;
//	rows = a.size();
//	cols = b.size();
//	printf("Rows = %d Cols = %d\n",rows,cols);
//	double overlap;
//#ifdef USE_VNL_HUNGARIAN
//	vnl_matrix<double> mat(rows,cols);
//	printf("Allocated Rows = %d Cols = %d\n",mat.rows(),mat.cols());
//	//int pa =0; 
//	//int pb =0;
//	for(int cr = 0; cr<rows; cr++)
//	{
//
//		bool enforce_overlap = false;
//		for(int cc=0; cc<cols; cc++)
//		{
//			overlap = features_box_overlap(a[cr],b[cc]);
//			if(overlap>1) // volume of overlap > 1
//			{
//				enforce_overlap = true;
//				break;
//			}
//		}
//		/*	if(rank==DEBUG_RANK)
//		{
//		if(enforce_overlap)
//		printf("%d/%d: I did enforce overlap for %d\n",rank,npes,cr);
//		else
//		printf("%d/%d: I did not enforce overlap for %d\n",rank,npes,cr);
//		}*/
//		for(int cc =0; cc<cols; cc++)
//		{
//			mat(cr,cc) = features_diff(a[cr],b[cc],enforce_overlap);
//		}
//	}
//	printf("About to call vnl_hungarian_algorithm\n");
//	vcl_vector<unsigned int> ret = vnl_hungarian_algorithm<double>(mat);
//	printf("Returned from vnl_hungarian_algorithm\n");
//	for(unsigned int counter=0; counter< ret.size(); counter++)
//	{
//		if(mxIsFinite(ret[counter]))
//		{
//			if(!mxIsFinite(mat(counter,ret[counter])))
//			{
//				ret[counter] = static_cast<unsigned int>(-1);
//			}
//		}
//	}
//	printf("Returning from getTimeAssociations\n");
//	return ret;
//#else
//	double * assignment  = (double*) malloc(rows*sizeof(double));
//	double * cost = (double*) malloc(sizeof(double));
//	double * distMatrixIn = (double*) malloc(rows *cols*sizeof(double));
//
//	if(assignment == NULL)
//		printf("Couldn't allocate memory assignment\n");
//
//	if(distMatrixIn == NULL)
//		printf("Couldn't allocate memory for distMatrixIn\n");
//
//	if(cost == NULL)
//		printf("Couldn't allocate memory for cost\n");
//
//	printf("%d/%d About to assign datamatrix values\n",rank,npes);
//	int pa =0; 
//	int pb =0;
//	for(int cr = 0; cr<rows; cr++)
//	{
//
//		//if(rank==0)
//		//{
//		//	printf("0/124: bbox %d %d %d %d %d %d \n",a[cr].bbox.sx,a[cr].bbox.sy,a[cr].bbox.sz,a[cr].bbox.ex,a[cr].bbox.ey,a[cr].bbox.ez);
//		//}
//		bool enforce_overlap = false;
//		for(int cc=0; cc<cols; cc++)
//		{
//			overlap = features_box_overlap(a[cr],b[cc]);
//			if(overlap>1) // volume of overlap > 1
//			{
//				enforce_overlap = true;
//				break;
//			}
//		}
//		//if(rank==DEBUG_RANK)
//		{
//			if(enforce_overlap)
//				printf("%d/%d: I did enforce overlap for %d\n",rank,npes,cr);
//			else
//				printf("%d/%d: I did not enforce overlap for %d\n",rank,npes,cr);
//		}
//		for(int cc =0; cc<cols; cc++)
//		{
//			distMatrixIn[cr+rows*cc] = features_diff(a[cr],b[cc],enforce_overlap);
//			//if(rank==DEBUG_RANK)
//			{
//				//	if(cr==0)
//				{
//					printf("distmatrix[%d,%d]=%0.3lf\n",cr,cc,distMatrixIn[cr+rows*cc]);
//				}
//			}
//		}
//
//	}
//	printf("%d/%d About to call assignmentsuboptimal1\n",rank,npes);
//	assignmentsuboptimal1(assignment,cost,distMatrixIn,rows,cols);
//	printf("%d/%d Exited assignmentsuboptimal1\n",rank,npes);
//	vcl_vector<unsigned int> vec(rows);
//	int assigned_count=0;
//	//	if(rank==DEBUG_RANK)
//	//	{
//	//		for(int counter=0; counter< rows; counter++)
//	//		{
//	//			printf("%d/%d assignment[%d] = %0.3lf\n",rank,npes,counter,assignment[counter]);
//	//		}
//	//	}
//	printf("%d/%d About to start assigning vec values\n",rank,npes);
//	for(int cr = 0; cr<rows; cr++)
//	{
//		if(assignment[cr]>-0.1)
//		{
//			vec[cr]=static_cast<unsigned int>(assignment[cr]+0.5);
//			//			if(rank==DEBUG_RANK)
//			//				printf("%d/%d : assigned_nums[%d] = %d\n",rank,npes,cr,vec[cr]);
//			assigned_count++;
//		}
//		else
//		{
//			vec[cr]=static_cast<unsigned int>(-1);
//		}
//	}
//	//	if(rank==DEBUG_RANK)
//	//	{
//	//		printf("%d/%d: assigned_count = %d\n",rank,npes,assigned_count);
//	//	}
//	//free(assignment); FIXME : was giving glibc corruption errors;
//	//free(cost);
//	//free(distMatrixIn);
//	return vec;
//
//#endif
//
//}

void getArrayFromStdVector(std::vector<FeaturesType> &f, FeaturesType	*&farray)
{
	farray = new FeaturesType[f.size()];
	for(unsigned int counter=0; counter<f.size(); counter++)
	{
		farray[counter]=f[counter];
	}
}

void getStdVectorFromArray(FeaturesType *farray, int n,std::vector<FeaturesType> &f)
{
	f.reserve(n);
	for(int counter=0; counter<n; counter++)
	{
		f.push_back(farray[counter]);
	}
}




//#define CACHE_PREFIX "D:/ucb dataset/output/ena/cache"
#define CACHE_PREFIX "cache"
//#define OUTPUT_PREFIX "output"
#define MAX_TIME 50
#define MAX_TAGS 4
#define MAX_LABEL 10000
#define VESSEL_CHANNEL 4 // FIXME : make it dynamic based on user input
#define PAUSE {printf("%d:>",__LINE__);scanf("%*d");}

bool compare(FeaturesType a, FeaturesType b)
{
	return a.time<b.time;
}
void createTrackFeatures(std::vector<FeaturesType> fvector[MAX_TIME][MAX_TAGS], std::vector<ftk::TrackFeatures> &tfs, int c,int num_t)
{
	int max_track_num = 0;
	for(int t = 0; t< num_t; t++)
	{
		for(unsigned int counter=0; counter< fvector[t][c-1].size(); counter++)
		{
			max_track_num = MAX(max_track_num,fvector[t][c-1][counter].num);
		}
	}

	for(int counter=1; counter <= max_track_num; counter++)
	{
		ftk::TrackFeatures trackf;
		trackf.intrinsic_features.clear();
		for(int t = 0; t< num_t;t++)
		{
		for(unsigned int counter1 = 0; counter1 < fvector[t][c-1].size(); counter1++)
		{
			if(fvector[t][c-1][counter1].num == counter)
			{
				trackf.intrinsic_features.push_back(fvector[t][c-1][counter1]);
			}
		}
		}
		std::sort(trackf.intrinsic_features.begin(),trackf.intrinsic_features.end(),compare);
		tfs.push_back(trackf);
		//PRINTF("Added %d elements to tfs\n",counter);
	}
}
int main(int argc, char **argv)
{
//	printf("Entered Main\n");
	
	/*int num_files = atoi(argv[1]);
	int num_channels= atoi(argv[2]);
	int counter = atoi(argv[3]);
	FILE * fp = fopen(argv[4],"r");
*/
	/*printf("counter =%d\n",counter);
	printf("channles =%d\n",num_channels);
	printf("files =%d\n",num_files);
	printf("filenames = %s\n", argv[4]);
	*/

	InputImageType::Pointer im[50];
	InputImageType::Pointer om[50];
		 

	/*for(int i=0;i<(2*num_files);i++)
	{
		fgets(file_name[i],1024,fp);
		file_name[i][strlen(file_name[i])-1]= '\0';
		printf("files:%s\n",file_name[i]);
		if(feof(fp))
			break;
	}*/

  
	//sprintf(unmix_cache, "%s/%s",CACHE_PREFIX,file_name[counter+num_files]);
  bool fexists = true;
  int num_channels = (argc-1)/2;
//  printf("I found argc = %d\n",argc);
 
  for(int counter=num_channels+1; counter < argc; counter++)
    if(!file_exists(argv[counter]))
      {
        fexists = false;
        break;
      }
	if(!fexists)
	{
		for(int c=1;c<=num_channels; c++)
		{
			im[c-1] = readImage<InputImageType>(argv[c]);
		}

		unmix_median(im,om,num_channels);
		for(int c=1;c<=num_channels; c++)
		{
			im[c-1]=om[c-1];
		}

		for(int c=1; c<= num_channels; c++)
		{
			writeImage<InputImageType>(im[c-1],argv[c+num_channels]);
		}
	}
	printf("End of Main");

	return 0;
}
