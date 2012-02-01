/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

//#include <stdio.h>
#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>
#include "find_median.cpp"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkDanielssonDistanceMapImageFilter.h>


using namespace std;

//#define COMPUTE_L2

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

// The following lines help us access the linear arrays raster and p like a multidimensional array
#define MAT(a,b,c) raster[(a)*npixels+(b)*w+(c)]
#define P(a,b,c) p[(a)*pimsize+(b)*pwidth+(c)]
#define M(a,b,c) matrix[(a)*pimsize+(b)*pwidth+(c)]
#define C(a,b,c) checked[(a)*pimsize+(b)*pwidth+(c)]
//#ifdef COMPUTE_L2
#define L(a,b,c) lmatrix[(a)*pimsize+(b)*pwidth+(c)]
//#endif

// window - max window width for X,Y directions
int window=2;
// window1 - max window width for Z direction
int window1=2;
// lower threshold
double lambda_lower_threshold = 18;//original 5,  changed  on 01,29,2012
// higher threshold
double lambda_higher_threshold = 23;//original 8,  changed  on 01,29,2012
//prune value
double prune=3;
double epsilonw=1.5,alpha1;
#define CFH(a,b) (((unsigned int)a)<<(b*4))

int psize,pimsize,pwidth ,pheight,pdepth ;
double planes[42][3];
unsigned char * null,*fore,*back;
unsigned char *p;
#define DEBUG_ printf
//compare function for sorting the data
int compare(const void *a,const void *b)
{
	return (int)(*(unsigned char*)a)-(int)(*(unsigned char*)b);
}
int compare_double(const void *a,const void *b)
{
	return (*(double*)a - *(double*)b);
}
struct data{
	int x,y,z;
};

struct return_data{
	int M;
	double lvalue;
};
void ParseArguments(int argc, char **argv)
{
	for (int counter=3; counter<argc; counter++)
	{
		istringstream s(argv[counter]);
		//char name[100];
		char ch;
		int d;
		sscanf(argv[counter],"%c=%d",&ch,&d);
		switch(ch)
		{
		case 'x': case 'y': window = d;
			DEBUG_("window = %d\n",d);
			break;
		case 'z':	window1 = d;
			DEBUG_("window1 = %d\n",d);
			break;
		case 'l':	lambda_lower_threshold =d;
			DEBUG_("lambda_lower_threshold = %d\n",d);
			break;
		case 'h':	lambda_higher_threshold=d;
			DEBUG_("lambda_higher_threshold = %d\n",d);
			break;
		case 'p':	prune=d;
			DEBUG_("prune = %d\n",d);
			break;
		case 'w': epsilonw=d/2.0;
			DEBUG_("epsilonw = %lf\n",d/2.0);
			break;
		default:
			printf("Parse error in %s\n",argv[counter]);
			return;
		}
	}
	printf("\n\nSuccessfully parsed %d parameters\n\n\n",argc-1);
}

void set_limits(int &zl, int &zu, int &yl, int &yu, int &xl, int&xu,int argc, char**argv)
{
	zl = atoi(argv[argc-6]);
	zu = atoi(argv[argc-5]);
	yl = atoi(argv[argc-4]);
	yu = atoi(argv[argc-3]);
	xl = atoi(argv[argc-2]);
	xu = atoi(argv[argc-1]);

}
int main(int argc, char**argv)
{
	time_t t1;
	t1 = time(NULL);
	unsigned char * raster;
	//unsigned char *maxintensity;
	unsigned int rwidth,rlength,rdepth;
	// pic file name

	printf("%d\n",argc);
	//if(argc >2)
	ParseArguments(argc,argv);
	//char picfile[1024];// dont segfault me plz!
	//printf("%s\n",argv[1]);
	char filenamebuff[1024];//="vessel.pic";
//	sprintf(filenamebuff,"%s.pic",argv[1]);
	typedef itk::Image<unsigned char,3> ImageType;
	typedef itk::ImageFileReader<ImageType> FileReaderType;

	FileReaderType::Pointer reader = FileReaderType::New();
	reader->SetFileName(argv[1]);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 1;
	}

	ImageType::Pointer im = reader->GetOutput();

	std::cout << "Image Loaded";

	rwidth = im->GetLargestPossibleRegion().GetSize()[0];
	rlength = im->GetLargestPossibleRegion().GetSize()[1];
	rdepth = im->GetLargestPossibleRegion().GetSize()[2];

	ImageType::Pointer imout = ImageType::New();
	imout->SetRegions(im->GetLargestPossibleRegion());
	imout->Allocate();

	/*
	FILE * fpi = fopen(filenamebuff,"rb");
	if(fpi == NULL)
	{
		printf("I tried this file '%s'\n",filenamebuff);
		printf("Couldn't open the file. Please check the filename and try again\n");
		return 0;
	}

	unsigned char *temp=new unsigned char[77];
	fread(temp,sizeof(unsigned char),76,fpi);

	rwidth = CFH(temp[0],0)+CFH(temp[1],2);
	rlength = CFH(temp[2],0)+CFH(temp[3],2);
	rdepth = CFH(temp[4],0)+CFH(temp[5],2);
	*/

	printf("I found rdepth, rlength, rwidth %d %d %d\n",rdepth,rlength,rwidth);
	//scanf("%*d");
	int npixels = rwidth*rlength;
	int w = rwidth;
	raster = (unsigned char*)malloc(npixels*rdepth*sizeof(unsigned char));
	if(raster==NULL)
	{
		printf("memory problem: couldn't allocate enough memory\n");
		_exit(0);
	}
	memcpy(raster,im->GetBufferPointer(),npixels*rdepth*sizeof(unsigned char));


	char buff [1024];// = "E:\\Arun\\My matlab codes\\deconvolved\\deconv";
	strcpy(buff,filenamebuff);


	const int skip  = 1;
	char filenameoutbuff[1024];
	sprintf(filenameoutbuff,"%s.npts",argv[1]);
	FILE * fp = fopen(filenameoutbuff,"w");
	
	double *lmatrix;
	int * matrix;
	bool *checked;
	double *l1matrix;
	double *l2matrix;
	const int gridsize_x = 1;
	const int gridsize_y = 1;
	const int gridx = 0;
	const int gridy = 0;
	alpha1 = 0.3;

	vector <data> queuep;
	queuep.reserve(200000);
	int neighbors[26][3] = {
		{-1,    -1,    -1},
		{	-1,     0,    -1},
		{	-1,     1,    -1},
		{	0,    -1,    -1},
		{	0,     0,    -1},
		{	0,     1,    -1},
		{	1,    -1,    -1},
		{	1,     0,    -1},
		{	1,     1,    -1},
		{	-1,    -1,     0},
		{	-1,     0,     0},
		{	-1,     1,     0},
		{	0,    -1,     0},
		{	0,     1,     0},
		{	1,    -1,     0},
		{	1,     0,     0},
		{	1,     1,     0},
		{	-1,    -1,     1},
		{	-1,     0,     1},
		{	-1,     1,     1},
		{	0,    -1,     1},
		{	0,     0,     1},
		{	0,     1,     1},
		{	1,    -1,     1},
		{	1,     0,     1},
		{	1,     1,     1}
	};

	int zl = 0,zu=rdepth-1,yl=gridy*(rlength/gridsize_y),yu=(gridy+1)*(rlength/gridsize_y)-1,xl=gridx*(rwidth/gridsize_x),xu=(gridx+1)*(rwidth/gridsize_x)-1;

	printf("zl zu yl yu xl xu = %d %d %d %d %d %d\n",zl,zu,yl,yu,xl,xu);

	// We find the various dimensions of the image for our ease of use
	psize = (zu-zl+1+2*window1)*(yu-yl+1+2*window)*(xu-xl+1+2*window);
	pimsize = (yu-yl+1+2*window)*(xu-xl+1+2*window);
	pwidth = (xu-xl+1+2*window);
	pheight = (yu-yl+1+2*window);
	pdepth = (zu-zl+1+2*window1);
	int gamma_size = (2*window+1)*(2*window1+1)*(2*window+1);
	// allocate memory and copy the data into a new array
	p = (unsigned char *)malloc(psize*sizeof(unsigned char));
	matrix = (int *) malloc(psize*sizeof(int));
	memset(matrix,0,psize*sizeof(int));
	printf("did I cross that?\n");
	checked = (bool *) malloc(psize*sizeof(bool));
	l1matrix = (double *)malloc(gamma_size*sizeof(double));
	l2matrix = (double *)malloc(gamma_size*sizeof(double));
	//#ifdef COMPUTE_L2
	lmatrix = (double *)malloc(psize*sizeof(double));
	if(l1matrix == NULL || l2matrix == NULL || lmatrix == NULL || matrix == NULL)
	{
		printf("mem alloc problem %p %p %p %p\n",l1matrix,l2matrix,lmatrix,matrix);

	}
	else
	{
		printf("No mem alloc problem\n");
	}
	//#endif
	for(int coz = 0; coz <pdepth; coz++)
		for(int coy = 0 ; coy <pheight; coy++)
			for(int cox = 0; cox <pwidth; cox++)
				M(coz,coy,cox)=0;

	printf("Problem here? pwidth = %d pheight = %d pdepth = %d psize = %d \n",pwidth, pheight, pdepth,psize);
	// replicated
	int count = 0;
	for (int coz=window1;coz<pdepth-window1;coz++)
		for(int coy=window;coy<pheight-window;coy++)
			for(int cox=window;cox<pwidth-window;cox++)
			{
				if(MAT(zl+coz-window1,yl+coy-window,xl+cox-window)>0)
				{
					count++;
				}
				P(coz,coy,cox)=MAT(zl+coz-window1,yl+coy-window,xl+cox-window);
			}

			printf("Problem here? pwidth = %d pheight = %d pdepth = %d\n",pwidth, pheight, pdepth);
			//'replicate' the boundaries
			for( int counter = 0; counter < window; counter++)
			{
				for(int counter1 =0; counter1 < pdepth; counter1++)
					for(int counter2 = 0; counter2 < pwidth; counter2++)
					{
						P(counter1,counter,counter2)=P(counter1,window,counter2);
						P(counter1,pheight-counter-1,counter2)=P(counter1,pheight-window-1,counter2);
					}
					for(int counter1 =0; counter1 < pdepth; counter1++)
						for(int counter2 = 0; counter2 < pheight; counter2++)
						{
							P(counter1,counter2,counter)=P(counter1,counter2,window);
							P(counter1,counter2,pwidth-counter-1)=P(counter1,counter2,pwidth-window-1);
						}
			}
			for( int counter = 0; counter < window1; counter++)
			{
				for(int counter1 =0; counter1 < pheight; counter1++)
					for(int counter2 = 0; counter2 < pwidth; counter2++)
					{
						P(counter,counter1,counter2)=P(window1,counter1,counter2);
						P(pdepth-counter-1,counter1,counter2)=P(pdepth-window-1,counter1,counter2);
					}
			}
			printf("mat >0 in %d points\n",count);
			printf("Did I come till epsilon?\n");
			// continue adding more variables
			//						double alpha1 = 0.3;
			//double epsilon = 1e-2;

			//double scaling_matrix[] = {window1, window, window};

			double planes_21[21][3]= {
        {1, 0, -1},
				{1, 0, 1},
				{1, -1, 0},
				{1, 1, 0},
				{0, 1, 1},
				{0, -1, 1},
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1},
				{0.5, 0.5, 1},
				{0.5, -0.5, 1},
				{-0.5, -0.5, 1},
				{-0.5, 0.5, 1},
				{1, 0.5, 0.5},
				{1, 0.5, -0.5},
				{1, -0.5, 0.5},
				{1, -0.5, -0.5},
				{0.5, 1, 0.5},
				{-0.5, 1, 0.5},
				{0.5, 1, -0.5},
				{-0.5, 1, -0.5}
        };

			// Reflect and store in a new array
			for(int counter = 0; counter < 21; counter++)
			{
				planes[counter][0]=planes_21[counter][0];
				planes[counter][1]=planes_21[counter][1];
				planes[counter][2]=planes_21[counter][2];
				planes[counter+21][0]=-planes_21[counter][0];
				planes[counter+21][1]=-planes_21[counter][1];
				planes[counter+21][2]=-planes_21[counter][2];
			}
			for(int counter =0; counter < 21; counter++)
			{
				double sum  = sqrt(planes[counter][0]*planes[counter][0]+planes[counter][1]*planes[counter][1]+planes[counter][2]*planes[counter][2]);
				planes[counter][0]/=sum;
				planes[counter][1]/=sum;
				planes[counter][2]/=sum;
			}
			//unsigned int gamma_size = (2*window1+1)*(2*window+1)*(2*window+1);
			double t[3]; // stores the lambda values

			null = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
			fore = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
			back = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
//			like = (double *)malloc((2*window1+1)*(2*window+1)*(2*window+1)*sizeof(float));
			//epsilonw=2.5;
			printf("epsilongw = %lf\n",epsilonw);
			queuep.clear();
			double best_t1,best_t0;
			double likelihood=0;
			for(int coz = window1;coz<pdepth-window1;coz+=skip)
			{
				for(int coy = window; coy<pheight-window;coy+=skip)// skip is used for sampling. We dont sample in Z direction
					for(int cox = window; cox<pwidth-window;cox+=skip)
					{
						best_t1 = best_t0 = -1;
						//the main for loop
						double L1 = -5000000;
						double L;
						double amax,bmax,cmax;
						//int pcount =0;
						int pccn =0;

						for(int wz=-window1;wz<=window1;wz++)
							for(int wy=-window;wy<=window;wy++)
								for(int wx=-window;wx<=window;wx++)
								{
									null[pccn++]=P(coz+wz,coy+wy,cox+wx);
								}

						//defined in stdlib to sort sequential arrays
						//	qsort(null,pccn,sizeof(unsigned char),compare);
						// find the median of the null hypothesis
						/*double sum_sort =0;
						int sum_sort_num =0;
						for(int counter_s=alpha1*pccn; counter_s<pccn-alpha1*pccn;counter_s++,sum_sort_num++)
						sum_sort+=null[counter_s];
						t[2]=sum_sort/sum_sort_num;
						*/
						//if(pccn%2==0)
						//	t[2]=(null[pccn/2]+null[pccn/2-1])/2.0;
						//else
						//	t[2]=null[pccn/2];

						t[2] = median(pccn,null);

						// a small threshold fixed to accelerate the program. Anything with a median of less than 2 is ignored for DEBUG

						if(t[2]<prune)
							continue;
						double a,b,c;
						//printf("I did come here\n");
						// now for the alternate hypothesis
						for(int co = 0; co < 21;co++)
						{

							a = planes[co][0];
							b = planes[co][1];
							c = planes[co][2];
							int pccf =0;
							int pccb =0;
							for(int wz =-window1;wz<=window1;wz++)
								for(int wy=-window;wy<=window;wy++)
									for(int wx=-window;wx<=window;wx++)
									{
										//FIXME
										if((wz*a+wy*b+wx*c)<= epsilonw  && (wz*a+wy*b+wx*c)>=-epsilonw)
											fore[pccf++]=P(coz+wz,coy+wy,cox+wx);
										else
											back[pccb++]=P(coz+wz,coy+wy,cox+wx);
									}

									//this is an inefficient way to find the median. 
									//I sort the arrays and find the middle element
									//qsort(fore,pccf,sizeof(unsigned char),compare);
									//qsort(back,pccb,sizeof(unsigned char),compare);
									//find the median of back
									//if(pccb%2==0)
									//	t[0]=(back[pccb/2]+back[pccb/2-1])/2.0;
									//else
									//	t[0]=back[pccb/2];
									////find the median of fore
									//if(pccf%2==0)
									//	t[1]=(fore[pccf/2]+fore[pccf/2-1])/2.0;
									//else
									//	t[1]=fore[pccf/2];

									t[0] = median(pccb,back);
									t[1] = median(pccf,fore);
								
									likelihood=0;
									if(t[1]-t[0]>=lambda_lower_threshold)
									{
										L=t[1]-t[0];

										if(t[2]==0)
										{
											likelihood = 1000000;//set some high value because likelihood is \inf
										}
										else
										{
											//new 
											likelihood += (t[1]+t[0]-2*t[2])/t[2]+log(t[2]/t[1]);
											if(t[0]!=0)
											{
												likelihood +=log(t[2]/t[0]);
											}
											if(likelihood<0)
											{
												likelihood = 1000000;
											}
										}
									}
									else
										likelihood=-5000000;
									if(likelihood>L1)
									{
										amax = a;
										bmax = b;
										cmax = c;
										best_t1 = t[1];
										best_t0 = t[0];
										L1 = likelihood;
										//if(L1 > 0)
										//	break;
									}
						}
						C(coz,coy,cox)=false;
						if(L1>0) // tau = 0
						{
							//printf("I came here\n");
							M(coz,coy,cox) = int(best_t1-best_t0+0.5);
							if(best_t1-best_t0 >= lambda_higher_threshold)
							{
								C(coz,coy,cox)=true;
								data d ;
								d.x = cox;
								d.y = coy;
								d.z = coz;
								queuep.push_back(d);
							}
							L(coz,coy,cox)=L1;
							//compute L2 for the point
							//double L2 = 0;

						}
						else
						{
							M(coz,coy,cox)=-1;


						}
					}
					printf("%s: %d coz\n",filenamebuff,coz);
					printf("%s: %ld seconds have elapsed since starting\n",filenamebuff, time(NULL)-t1);
					fflush(stdout);
			}
			// now start with seeds and expand at all places again
			// Things TODO :
			// 
			//free(p);
			int pc = -1;
			
			while(1)
			{
				pc++;
				if((unsigned int)pc == queuep.size())
					break;
				data seedp = queuep[pc];
				assert(M(seedp.z,seedp.y,seedp.x)>0);
				int l1 = M(seedp.z,seedp.y,seedp.x);

				//double l2 = L(seedp.z,seedp.y,seedp.x);

				for (int co = 0; co<26; co++) // 3x3x3 neighbors
				{
					int x1 = seedp.x+neighbors[co][0];
					int y1 = seedp.y+neighbors[co][1];
					int z1 = seedp.z+neighbors[co][2];
					if(!(x1>=window && y1>=window && z1>=window1 && x1<pwidth-window && y1<pheight-window && z1 <pdepth-window1))
					{
						continue;
					}
					if(C(z1,y1,x1)==1)
						continue;
					
					if(M(z1,y1,x1)>=MAX(l1-8,lambda_lower_threshold))
					{
						data D;
						D.x = x1;
						D.y = y1;
						D.z = z1;
						queuep.push_back(D);
						C(z1,y1,x1)=1;
					}
				}
			}
			printf("About to print the data\n");
			ImageType::IndexType output_index;
			for(int coz = window1;coz<pdepth-window1;coz++)
				for(int coy = window; coy<pheight-window;coy+=1)// skip is used for sampling sometimes. We dont sample in Z direction
					for(int cox = window; cox<pwidth-window;cox+=1)
					{
						if(C(coz,coy,cox)==1)
						{
							//int minimum_x=MAX(1,cox-window),minimum_y=MAX(1,coy-window),minimum_z=MAX(1,coz-window1);
							//int maximum_x=MIN(rwidth/gridsize_x,cox+1-window),maximum_y=MIN(rlength/gridsize_y,coy+1-window),maximum_z=MIN(rdepth,coz+1-window1);

							for(int cx = cox+1-window;cx<=cox+1-window;cx++)
								for(int cy = coy+1-window;cy<=coy+1-window;cy++)
								{
									output_index[0] = xl+cx-1;
									output_index[1] = yl+cy-1;
									output_index[2] = zl+coz-window1;

									fprintf(fp,"%d %d %d %d",zl+coz+1-window1,yl+cy,xl+cx,M(coz,coy,cox));
									imout->SetPixel(output_index,255);
									//#ifdef COMPUTE_L2
									fprintf(fp," %lf",L(coz,coy,cox));
									//#endif
									fprintf(fp,"\n");
								}

						}
					}
			
			typedef itk::ImageFileWriter<ImageType> FileWriterType;
			FileWriterType::Pointer writer = FileWriterType::New();
			writer->SetFileName(argv[2]);
			writer->SetInput(imout);
			try
			{
				writer->Update();
			}
			catch (itk::ExceptionObject &e)
			{
				std::cout << e << std::endl;
				return 1;
			}

			char distance_map_file[1024];
			strcpy(distance_map_file,argv[2]);
			strcpy(&distance_map_file[strlen(distance_map_file)-4],"_distance_map.tif");
			char vector_map[1024];
			strcpy(vector_map,argv[2]);
			strcpy(&vector_map[strlen(vector_map)-4],"_vector_map.mhd");

			typedef itk::Image<short int,3> DistanceImageType;
			typedef itk::DanielssonDistanceMapImageFilter<ImageType,DistanceImageType> DistanceMapFilterType;
			typedef DistanceMapFilterType::VectorImageType OffsetImageType;

			DistanceMapFilterType::Pointer distfilter = DistanceMapFilterType::New();
			distfilter->SetInput(imout);
			distfilter->InputIsBinaryOn();
			try
			{
				distfilter->Update();
			}
			catch (itk::ExceptionObject &e)
			{
				std::cout << e << std::endl;
				return 1;
			}

			typedef itk::ImageFileWriter<DistanceImageType> DistFileWriter;
			DistFileWriter::Pointer dwrite = DistFileWriter::New();
			dwrite->SetFileName(distance_map_file);
			dwrite->SetInput(distfilter->GetOutput());
			try
			{
				dwrite->Update();
			}
			catch (itk::ExceptionObject &e)
			{
				std::cout << e << std::endl;
				return 1;
			}
			
			typedef itk::ImageFileWriter<OffsetImageType> OffsetFileWriter;
			OffsetFileWriter::Pointer offwriter = OffsetFileWriter::New();
			offwriter->SetFileName(vector_map);
			offwriter->SetInput(distfilter->GetVectorDistanceMap());
			try
			{
				offwriter->Update();
			}
			catch (itk::ExceptionObject &e)
			{
				std::cout << e << std::endl;
				return 1;
			}

			free(matrix);
			free(lmatrix);
			free(checked);
			free(l1matrix);
			free(l2matrix);
			free(raster);
			fclose(fp);
			printf("%ld seconds \n",time(NULL)-t1);
			return 0;
}
