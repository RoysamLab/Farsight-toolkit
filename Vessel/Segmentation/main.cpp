///*=========================================================================
//Copyright 2009 Rensselaer Polytechnic Institute
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License. 
//=========================================================================*/

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
#include "omp.h"
#include "find_median.cpp"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkDanielssonDistanceMapImageFilter.h"

//lets try tensor 66
#include "rtvl/rtvl_tensor.hxx"
#include "rtvl/rtvl_vote.hxx"
#include "rtvl/rtvl_votee.hxx"
#include "rtvl/rtvl_voter.hxx"
#include "rtvl/rtvl_weight_smooth.hxx"

//#include <rtvl_tensor.hxx>
//#include <rtvl_vote.hxx>
//#include <rtvl_votee.hxx>
//#include <rtvl_voter.hxx>
//#include <rtvl_weight_smooth.hxx>


#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"

#include "vcl_iostream.h"

//for vesselness
#include "itkCastImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkDivideImageFilter.h"

typedef itk::Image<unsigned char,3> ImageType;
typedef itk::Image<float,3> ImageTypeFloat;



using namespace std;

//#define COMPUTE_L2

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

// The following lines help us access the linear arrays raster and p like a multidimensional array
#define MAT(a,b,c) raster[(a)*npixels+(b)*w+(c)]
#define P(a,b,c) p[(a)*pimsize+(b)*pwidth+(c)]
#define M(a,b,c) matrix[(a)*pimsize+(b)*pwidth+(c)]
#define C(a,b,c) checked[(a)*pimsize+(b)*pwidth+(c)]
#define SM(a,b,c) sparmark[(a)*pimsize+(b)*pwidth+(c)]
#define Mvm1(a,b,c) Vm1[(a)*pimsize+(b)*pwidth+(c)]
#define Mvm2(a,b,c) Vm2[(a)*pimsize+(b)*pwidth+(c)]
#define CVOTEE(a,b,c) Cvotee[(a)*pimsize+(b)*pwidth+(c)]
//#define Que0(a,b,c) queuep0.at((a)*pimsize+(b)*pwidth+(c))
//#define MItin(a,b,c)   mark_iterin[(a)*pimsize+(b)*pwidth+(c)]
//#define MItout(a,b,c)  mark_iterout[(a)*pimsize+(b)*pwidth+(c)]
//#ifdef COMPUTE_L2
#define L(a,b,c) lmatrix[(a)*pimsize+(b)*pwidth+(c)]
//#endif

// window - max window width for X,Y directions
int window=5,window_fix=5;
// window1 - max window width for Z direction
int window1=3;
// lower threshold
double lambda_lower_threshold=11,lambda_lower_threshold_fix = 11;
// higher threshold
double lambda_higher_threshold=14,lambda_higher_threshold_fix = 14;
//prune value
double prune=3,prune_adaptiveL=12,prune_adaptiveH=120;
double epsilonw=1.5,alpha1;
///threshold for analyse decomposing results
double thre = 0.8;//0.4//0.5//0.6//0.8 best//1.0//1.8//2.8//2.4//2.2//<1.6 lead to complete
//sparsize for sparse voting
int sparsize=5,sparsize1=2;///10_2,5_2is the same but much quicker
const int itert=8;//set iteration time
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

//struct hypo{
//	vector<double> t;
//	vector<unsigned char> null,back,fore;
//	double amax,bmax,cmax,a,b,c;
//	double L1,likelihood,best_t1,best_t0;
//};

struct vmatrix{
	vnl_matrix_fixed<double, 3, 3> voting_matrix;
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

vector<data> sparse_mark(int sparsize,int sparsize1,bool * mark_input,bool * mark_output,vector<data> queuep,vector<data> que_output)
{
	int poz,poy,pox;
	printf("Start Neighbor Region Mark for Sparse Voting\n");
	printf("the input queue size is %d\n",queuep.size());
	que_output.clear();
	//vector <data> que_output(queuep);
	//que_output=queuep;
	//que_output.assign(queuep.begin(),queuep.end());
	//printf("que_output size is %d\n",que_output.size());
	bool * sparmark;
	sparmark = (bool *) malloc(psize*sizeof(bool));
	if (sparmark==NULL)
		printf("Fatal problem: could not malloc enough memory for sparmark\n");
	memcpy(sparmark,mark_input,psize*sizeof(bool));
	int control_spar=queuep.size();
	//printf("que_output size is %d\n",que_output.size());



	for (int qcounter=0;qcounter<control_spar;qcounter++)
	{
		poz=queuep[qcounter].z;
		poy=queuep[qcounter].y;
		pox=queuep[qcounter].x;
		/*if ((poz<window1)||(poy<window)||(pox<window))
			continue;*/
		for (int cz=MAX(poz-sparsize1,window1);cz<MIN(poz+sparsize1,pdepth-window1);cz++)
			for (int cy=MAX(poy-sparsize,window);cy<MIN(poy+sparsize,pheight-window);cy++)
				for (int cx=MAX(pox-sparsize,window);cx<MIN(pox+sparsize,pwidth-window);cx++)
				{
					if (SM(cz,cy,cx)==0)
					{
						data po;
						po.z=cz;
						po.y=cy;
						po.x=cx;
						que_output.push_back(po);
						SM(cz,cy,cx)=1;
					}

				}

	}
	memcpy(mark_output,sparmark,psize*sizeof(bool));
	free(sparmark);
	printf("sparse mark complete\n");
	printf("que_output size is %d\n",que_output.size());
	return que_output;


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
	time_t time_cal[itert+1];
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

	//itk::ImageRegionIterator<ImageType> imit(im,im->GetRequestedRegion());
	//#pragma omp parallel for num_threads(14)
	//for (imit=imit.Begin();imit<=imit.End();++imit)
	//{
	//}



	////for vesselness
	///////////////////////////////////////
	//typedef itk::CastImageFilter< ImageType, ImageTypeFloat > CastFilterType;
	//CastFilterType::Pointer castFilter = CastFilterType::New();
	//castFilter->SetInput(im);
	//castFilter->Update();
	//ImageTypeFloat::Pointer floatout=castFilter->GetOutput();
	//
	//typedef itk::StatisticsImageFilter<ImageTypeFloat> StatisticsImageFilterType;
	//StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New ();
	//statisticsImageFilter->SetInput(floatout);
	//statisticsImageFilter->Update();
	//float maximus=statisticsImageFilter->GetMaximum();



	//float sigma_min = 1.0f; //0.5f;
	//float sigma_max = 10.0f; //4.0f;
	//int sigma_intervals=10;
	////int sigma_steps = 5;

	//float alpha = 0.5, beta = 0.5, gamma = 0.25; //5.0//0.025;
	//alpha=alpha*maximus;
	//beta=beta*maximus;
	//gamma=gamma*maximus;

	////int obj_dim = objectness_type; //1; //0: Blobness, 1: Vesselness, 2: Plateness

	//MultiScaleHessianFilterType::Pointer multi_scale_Hessian = MultiScaleHessianFilterType::New();
	//multi_scale_Hessian->SetInput(floatout);
	//multi_scale_Hessian->SetSigmaMin(sigma_min);
	//multi_scale_Hessian->SetSigmaMax(sigma_max);
	//multi_scale_Hessian->SetNumberOfSigmaSteps(sigma_intervals);

	////ObjectnessFilterType::Pointer objectness_filter = ObjectnessFilterType::New();
	//ObjectnessFilterType::Pointer objectness_filter = multi_scale_Hessian->GetHessianToMeasureFilter();
	//objectness_filter->SetScaleObjectnessMeasure(false);
	//objectness_filter->SetBrightObject(true);
	//objectness_filter->SetAlpha(alpha);
	//objectness_filter->SetBeta(beta);
	//objectness_filter->SetGamma(gamma);
	//objectness_filter->SetObjectDimension(1);/////);
	////std::cout << obj_measures.alpha << std::endl << obj_measures.beta << std::endl << obj_measures.gamma << std::endl;

	//multi_scale_Hessian->Update();

	//ImageTypeFloat::Pointer floatptr = multi_scale_Hessian->GetOutput();

	//StatisticsImageFilterType::Pointer statisticsImageFilter2 = StatisticsImageFilterType::New ();
	//statisticsImageFilter2->SetInput(floatptr);
	//statisticsImageFilter2->Update();
	//float maximus2=statisticsImageFilter2->GetMaximum();
	//
	////itk::Index<3> idx;
	////idx[0]=55;
	////idx[1]=68;
	////idx[2]=21;
	////float heihei=floatptr->GetPixel(idx);
	////std::cout<<heihei<<std::endl;
	typedef itk::DivideImageFilter <ImageTypeFloat, ImageTypeFloat, ImageTypeFloat > DivideImageFilterType;
	//DivideImageFilterType::Pointer divideImageFilter = DivideImageFilterType::New();
	//divideImageFilter->SetInput1(floatptr);
	//divideImageFilter->SetInput2(maximus2);
	//divideImageFilter->Update();



	//char vesselnessout[1024];
	//sprintf(vesselnessout,"%s_vesselness.mhd",argv[1]);
	//typedef itk::ImageFileWriter<ImageTypeFloat> FileWriterTypeFloat;
	//FileWriterTypeFloat::Pointer floatwriter = FileWriterTypeFloat::New();
	//floatwriter->SetFileName(vesselnessout);
	//floatwriter->SetInput(divideImageFilter->GetOutput());
	//try
	//{
	//	floatwriter->Update();
	//}
	//catch (itk::ExceptionObject &e)
	//{
	//	std::cout << e << std::endl;
	//	return 1;
	//}
	///////////////////////////////////////////////

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
	//vector <data> queuep0;
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
			//double t[3]; // stores the lambda values

			/*null = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
			fore = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
			back = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));*/
//			like = (double *)malloc((2*window1+1)*(2*window+1)*(2*window+1)*sizeof(float));
			//epsilonw=2.5;
			printf("epsilongw = %lf\n",epsilonw);
			queuep.clear();
			double best_t1,best_t0;
			//double likelihood=0;
			int NeighborVectorLength=(2*window+1)*(2*window+1)*(2*window1+1);
			//int StructVectorLength=(pheight-2*window)*(pwidth-2*window)*(pdepth-2*window1);
			//vector<hypo> hypotest;
			//hypotest.resize(StructVectorLength);
			//queuep0.resize(Que0Length);
			/*vector<unsigned char> null,fore,back;
			null.resize(NeighborVectorLength);
			fore.resize(NeighborVectorLength);
			back.resize(NeighborVectorLength);
			vector<double> t;
			t.resize(3);*/

            #pragma omp parallel for num_threads(14)
			//#pragma omp parallel for 16
			//#pragma omp parallel for num_threads(10) schedule(dynamic, 1)
			for(int coz = window1;coz<pdepth-window1;++coz)//coz+=skip
			{
				//#pragma omp parallel for num_threads(10) schedule(dynamic, 1)
				for(int coy = window; coy<pheight-window;++coy)//coy+=skip skip is used for sampling. We dont sample in Z direction
					for(int cox = window; cox<pwidth-window;++cox)//cox+=skip
					{
						//hypo hy;
						//double t[3];
						vector<double> t;
						t.resize(3);
						/*null = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
						fore = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
						back = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));*/
						vector<unsigned char> null,fore,back;
						null.resize(NeighborVectorLength);
						fore.resize(NeighborVectorLength);
						back.resize(NeighborVectorLength);
						/*unsigned char *null = 0;
						null = new unsigned char[];*/
						/*unsigned char *null=new unsigned char[847];
						unsigned char *back=new unsigned char[847];
						unsigned char *fore=new unsigned char[847];*/
						double best_t1,best_t0;
						best_t1 = best_t0 = -1;
						//double t[3];

						//the main for loop
						double L1 = -5000000,likelihood=0;
						//double L;
						double amax,bmax,cmax;
						//int pcount =0;
						int pccn =0;

						for(int wz=-window1;wz<=window1;wz++)
							for(int wy=-window;wy<=window;wy++)
								for(int wx=-window;wx<=window;wx++)
								{
									/*null.at(pccn)=P(coz+wz,coy+wy,cox+wx);
									pccn++;*/
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
							
						/*int lengthnull=null.size();
						unsigned char *nullarray;
						nullarray=new unsigned char[lengthnull];
						for (int ii=0;ii<lengthnull;++ii)
						{
							nullarray[ii]=null.at(ii);
						}*/
								t[2] = median(pccn,null);
						// a small threshold fixed to accelerate the program. Anything with a median of less than 2 is ignored for DEBUG

						if(t[2]<prune)
							continue;
						//if(t[2]<prune_adaptiveL)//7
						//{
						//	lambda_lower_threshold=lambda_lower_threshold_fix-5;
						//	lambda_higher_threshold=lambda_higher_threshold_fix-5;

						//}
						//if(prune_adaptiveH<t[2])//10
						//{
						//	lambda_lower_threshold=lambda_lower_threshold_fix+2;
						//	lambda_higher_threshold=lambda_higher_threshold_fix+2;
						//	window=window_fix-2;
						//	//int pccn =0;
						//	//for(int wz=-window1;wz<=window1;wz++)
						//	//	for(int wy=-window;wy<=window;wy++)
						//	//		for(int wx=-window;wx<=window;wx++)
						//	//		{
						//	//			null[pccn++]=P(coz+wz,coy+wy,cox+wx);
						//	//		}


						//}

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
										{
											/*fore.at(pccf)=P(coz+wz,coy+wy,cox+wx);
											pccn++;*/
											fore[pccf++]=P(coz+wz,coy+wy,cox+wx);
										}
										else
										{
											/*back.at(pccb)=P(coz+wz,coy+wy,cox+wx);
											pccb++;*/
											back[pccb++]=P(coz+wz,coy+wy,cox+wx);
										}
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

									/*int lengthback=back.size();
									unsigned char *backarray;
									backarray=new unsigned char[lengthback];
									for (int ii=0;ii<lengthback;++ii)
									{
										backarray[ii]=back.at(ii);
									}
									int lengthfore=fore.size();
									unsigned char *forearray;
									forearray=new unsigned char[lengthfore];
									for (int ii=0;ii<lengthfore;++ii)
									{
										forearray[ii]=fore.at(ii);
									}*/
									t[0] = median(pccb,back);
									t[1] = median(pccf,fore);
								
									likelihood=0;
									if(t[1]-t[0]>=lambda_lower_threshold)
									{
										//L=t[1]-t[0];

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
								//data d ;
								//d.x = cox;
								//d.y = coy;
								//d.z = coz;
								//#pragma omp critical
								//{
								//	//queuep.push_back(d);
								//	Que0(coz,coy,cox)=d;
								//}
							}
							L(coz,coy,cox)=L1;
							//compute L2 for the point
							//double L2 = 0;

						}
						else
						{
							M(coz,coy,cox)=-1;


						}
						//lambda_lower_threshold=lambda_lower_threshold_fix;
						//lambda_higher_threshold=lambda_higher_threshold_fix;
						//window=window_fix;
						/*delete null;
						delete fore;
						delete back;*/
						null.clear();
						fore.clear();
						back.clear();
						t.clear();
						/*free(null);
						free(back);
						free(fore);*/
					}
					printf("%s: %d coz\n",filenamebuff,coz);
					printf("%s: %ld seconds have elapsed since starting\n",filenamebuff, time(NULL)-t1);
					fflush(stdout);//need revise for omp
			}
			// now start with seeds and expand at all places again
			// Things TODO :
			// 
			//free(p);
			for(int coz = window1;coz<pdepth-window1;++coz)//coz+=skip
				for(int coy = window; coy<pheight-window;++coy)//coy+=skip skip is used for sampling. We dont sample in Z direction
					for(int cox = window; cox<pwidth-window;++cox)
					{
						if (C(coz,coy,cox)==true)
						{
							data d ;
						    d.x = cox;
							d.y = coy;
							d.z = coz;
							queuep.push_back(d);
						}
					}
			int pc = -1;
			
			lambda_lower_threshold=lambda_lower_threshold_fix;
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
			printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);

			printf("About to print the data after initial segmentation\n");
			ImageType::Pointer imout = ImageType::New();
			imout->SetRegions(im->GetLargestPossibleRegion());
			imout->Allocate();
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
			printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
			time_cal[0]=time(NULL)-t1;

			//-------------------------------------now for tensor voting----------------------------------------------------
			//
			///////////////////////////////////////////////////////////////////////
			
			for (int iter=0; iter<itert; iter++)
			{
				printf("%d-th iteration start...\n",iter);
				//define ball tensor for initialization
				vnl_matrix_fixed<double, 3, 3> voter_matrix;
				//vnl_matrix_fixed<double, 3, 3> votee_matrix(0.0);
				voter_matrix(0,0) = 1;
				voter_matrix(0,1) = 0;
				voter_matrix(0,2) = 0;
				voter_matrix(1,0) = 0;
				voter_matrix(1,1) = 1;
				voter_matrix(1,2) = 0;
				voter_matrix(2,0) = 0;
				voter_matrix(2,1) = 0;
				voter_matrix(2,2) = 1;
				//vnl_vector_fixed<double,3> voter_location(0.0);
				//vmatrix temp;
				int control1;
				int control2;
				//int StructVectorLength=(pheight-2*window)*(pwidth-2*window)*(pdepth-2*window1)+1;
				//control=(pwidth-2*window)*(pheight-2*window)*(pdepth-2*window1)-1;
				vector<vmatrix> Vm1;//store the votee_matrix result for token refinement,corresponding to each point in queuep
				vector<vmatrix> Vm2;
				//vector<vmatrix> Vm3;
				//Vm.reserve(2000000);
				Vm1.clear();
				
				/*std::vector< std::vector < vmatrix > > test;
				test.resize(control1);*/

				printf("Initialization for tensor voting complete\n");



				//Token Refinement
				printf("Start Token Refinement...\n");
				//control1=queuep.size();
				Vm1.resize(psize);
				printf("the queuep size is %d\n",control1);
				//#pragma omp parallel for 16
				#pragma omp parallel for num_threads(14)
				//#pragma omp parallel for num_threads(14) schedule(dynamic, 1)
				for(int coz = window1;coz<pdepth-window1;++coz)//coz+=skip
					for(int coy = window; coy<pheight-window;++coy)//coy+=skip skip is used for sampling. We dont sample in Z direction
						for(int cox = window; cox<pwidth-window;++cox)
				{
					if (C(coz,coy,cox)!=1)
						continue;
					vnl_matrix_fixed<double, 3, 3> votee_matrixaccu1(0.0);
					vnl_vector_fixed<double, 3> votee_location;
					votee_location[0]=(double)cox;
					votee_location[1]=(double)coy;
					votee_location[2]=(double)coz;
					for (int cz=MAX(coz-sparsize1,window1);cz<MIN(coz+sparsize1,pdepth-window1);++cz)
						for (int cy=MAX(coy-sparsize,window);cy<MIN(coy+sparsize,pheight-window);++cy)
							for (int cx=MAX(cox-sparsize,window);cx<MIN(cox+sparsize,pwidth-window);++cx)//maybe should use different sparsize for token refinement
							{
								if (C(cz,cy,cx)!=1)
									continue;
								//data voterp = queuep[voterc];
								vnl_vector_fixed<double,3> voter_location(0.0);
								voter_location[0]=(double)cx;
								voter_location[1]=(double)cy;
								voter_location[2]=(double)cz;
								vnl_matrix_fixed<double, 3, 3> votee_matrix(0.0);
								// Use "rtvl_tensor" to decompose the matrix.		
								rtvl_tensor<3> voter_tensor(voter_matrix);

								// Use "rtvl_voter" to encapsulate a token (location + input tensor).
								rtvl_voter<3> voter(voter_location, voter_tensor);

								// Use "rtvl_votee" to encapsulate a site (location + output tensor).
								rtvl_votee<3> votee(votee_location, votee_matrix);
								//vcl_cout << vcl_endl;

								// Choose a weight profile, initialized with spatial scale.
								rtvl_weight_smooth<3> tvw(1.0);

								// Compute one vote.
								rtvl_vote(voter, votee, tvw);

								votee_matrixaccu1+=votee_matrix;
								//printf("%d-th voter for Token Refinement\n",voterc);


							}
							vmatrix temp;
							temp.voting_matrix=votee_matrixaccu1;
							Mvm1(coz,coy,cox)=temp;
							//Vm1.at(voteec)=temp;
							//Vm1.push_back(temp);
							//test.at(voteec).push_back(temp);
							//printf("%d-th votee complete for Token Refinement\n",voteec);
				}
				printf("Token Refinement complete\n");
				printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);

				///Sparse Voting Mark

				//int sparsize=10,sparsize1=2;
				vector<data> que_output;
				bool * mark_sparout;
				mark_sparout = (bool *) malloc(psize*sizeof(bool));
				if (mark_sparout==NULL)
					printf("Fatal problem: could not malloc enough memory for mark_sparout\n");
				que_output=sparse_mark(sparsize,sparsize1,checked,mark_sparout,queuep,que_output);
				//printf("que_output size is %d\n",que_output.size());
				printf("Sparse Voting region marked after Token Refinement\n");
				printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);

				///Sparse Voting
				//int itert=1;
				//control2=Vm1.size();
				vector <data> result(que_output);
				//result=que_output;
				//vector <data> queue_voter(queuep);
				//vector <data> queue_votee(que_output);
				que_output.clear();

				bool * Cvotee;
				Cvotee = (bool *) malloc(psize*sizeof(bool));
				if (Cvotee==NULL)
				printf("Fatal problem: could not malloc enough memory for Cvotee\n");
				memcpy(Cvotee,mark_sparout,psize*sizeof(bool));
				free(mark_sparout);

				//bool * mark_iterout;
				//mark_iterout = (bool *) malloc(psize*sizeof(bool));
				//if (mark_iterout==NULL)
				//printf("Fatal problem: could not malloc enough memory for mark_iterout\n");

				printf("Start Sparse Voting...\n");
				//int voteenum;
				//vector<data> queueacc;
				//queueacc.clear();
				//queueacc=queuep;

				//for (int iter=0; iter<itert; iter++)
				//{
				//vector <data> queue_votee(que_output);
				//printf("%d-th iteration begins...\n",iter);
				//queue_votee=que_output;
				//voteenum=queue_votee.size();
				/*int lenGth;
				lenGth=queue_votee.size();
				for (int it=0;it<lenGth;it++)
				{
				data DAta=queue_votee[it];
				queueacc.push_back(DAta);

				}*/
				Vm2.clear();
				///////////main loop for votee, one voter to one votee at a time						
				//int vcounter=0;//store matrix for votee
				Vm2.resize(psize);
				//#pragma omp parallel for 16
				#pragma omp parallel for num_threads(14)
				//#pragma omp parallel for num_threads(14) schedule(dynamic, 1)
				for(int coz = window1;coz<pdepth-window1;++coz)//coz+=skip
					for(int coy = window; coy<pheight-window;++coy)//coy+=skip skip is used for sampling. We dont sample in Z direction
						for(int cox = window; cox<pwidth-window;++cox)
				{
					//printf("did i come here 66_1\n");
					if (CVOTEE(coz,coy,cox)!=1)
						continue;
					vnl_matrix_fixed<double, 3, 3> votee_matrixaccu2(0.0);
					//vote to one point each time
					vnl_vector_fixed<double, 3> votee_location;
					//int tempx,tempy,tempz;
					//tempx=queue_votee[ivotee].x;
					//tempy=queue_votee[ivotee].y;
					//tempz=queue_votee[ivotee].z;
					votee_location[0] = (double)cox;
					votee_location[1] = (double)coy;
					votee_location[2] = (double)coz;
					//printf("the value of control2 is %d\n",control2);

					/////loop for voter
					for (int cz=MAX(coz-sparsize1,window1);cz<MIN(coz+sparsize1,pdepth-window1);++cz)
						for (int cy=MAX(coy-sparsize,window);cy<MIN(coy+sparsize,pheight-window);++cy)
							for (int cx=MAX(cox-sparsize,window);cx<MIN(cox+sparsize,pwidth-window);++cx)
					{
						//printf("did i come here 66_2\n");
						if (C(cz,cy,cx)!=1)
							continue;
						vnl_vector_fixed<double,3> voter_location(0.0);
						/*int tx,ty,tz;
						tx=queue_voter[voterc].x;
						ty=queue_voter[voterc].y;
						tz=queue_voter[voterc].z;*/
						voter_location[0]=(double)cx;
						voter_location[1]=(double)cy;
						voter_location[2]=(double)cz;
						/*if (abs(tempx-tx)>sparsize)
							continue;
						if (abs(tempy-ty)>sparsize)
							continue;
						if (abs(tempz-tz)>sparsize1)
							continue;*/

						vnl_matrix_fixed<double, 3, 3> votee_matrix(0.0);


						//voter_matrix = Vm1[voterc].voting_matrix;
						voter_matrix = Mvm1(cz,cy,cx).voting_matrix;

						// Use "rtvl_tensor" to decompose the matrix.		
						rtvl_tensor<3> voter_tensor(voter_matrix);

						// Use "rtvl_voter" to encapsulate a token (location + input tensor).
						rtvl_voter<3> voter(voter_location, voter_tensor);

						// Use "rtvl_votee" to encapsulate a site (location + output tensor).
						rtvl_votee<3> votee(votee_location, votee_matrix);
						//vcl_cout << vcl_endl;

						// Choose a weight profile, initialized with spatial scale.
						rtvl_weight_smooth<3> tvw(1.0);

						// Compute one vote. 
						rtvl_vote(voter, votee, tvw);

						votee_matrixaccu2+=votee_matrix;
						//printf("%d-th voter for tensor",voterc);
						//printf("did i come here 66_4\n");

					}
					//printf("did i come here 66_3\n");
					vmatrix temp;
					temp.voting_matrix=votee_matrixaccu2;
					Mvm2(coz,coy,cox)=temp;
					//Vm2.at(ivotee)=temp;
					//Vm2.push_back(temp);
					//vcounter++;
					//printf("%d-th votee complete for Sparse Voting\n",vcounter);

				}
				free(Cvotee);
				//result=queue_votee;
				//printf("did i come here 66_4\n");
				//printf("%d-th iteration for tensor",iter);
				//if (iter==itert-1)
				//{
				//	//result=queueacc;
				//	result=queue_votee;
				//}
				//control2=Vm2.size();
				//vector <data> queue_voter(queue_votee);
				//queue_voter=queue_votee;
				//Vm1.clear();
				//Vm1=Vm2;
				//call sparse mark function to expand the region
				//que_output=sparse_mark(sparsize,sparsize1,mark_iterin,mark_iterout,queue_votee,que_output);
				//memcpy(mark_iterin,mark_iterout,psize*sizeof(bool));				
				//queue_votee.clear();

				//if (iter<itert-1)
				// vector <data> queue_votee(que_output);
				//vector <data> result(queue_votee);
				//printf("%d-th iteration complete\n",iter);
				//printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);

				//}
				//free(mark_iterin);
				//free(mark_iterout);
				Vm1.clear();
				printf("Sparse Voting complete\n");
				printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
				//Vm2.clear();

				// Decompose the result for sparse voting.
				printf("Start Decomposing the result for sparse voting\n");

				char saliencymap[1024];
				sprintf(saliencymap,"%s_saliencymap_%d_iteration.mhd",argv[2],iter);
				ImageTypeFloat::Pointer imsalout = ImageTypeFloat::New();
				imsalout->SetRegions(im->GetLargestPossibleRegion());
				imsalout->Allocate();
				ImageTypeFloat::IndexType imsalout_index;

				double maxsaliency=0;
				double sums=0;
				double means;
				double dividerv=0;
				double crit_0,crit_1,crit_2;
				int voteesize=result.size();
				int coz,coy,cox;
				vnl_matrix_fixed<double, 3, 3> votee_matrix(0.0);
				for (int voteecount=0;voteecount<voteesize;voteecount++)
				{
					coz=result[voteecount].z;
					coy=result[voteecount].y;
					cox=result[voteecount].x;
					if (C(coz,coy,cox)==1)
					{
						continue;
					}
					imsalout_index[0]=cox-window;
					imsalout_index[1]=coy-window;
					imsalout_index[2]=coz-window1;
					//votee_matrix=Vm2[voteecount].voting_matrix;
					votee_matrix=Mvm2(coz,coy,cox).voting_matrix;
					rtvl_tensor<3> votee_tensor(votee_matrix);
					crit_0=votee_tensor.saliency(0);
					imsalout->SetPixel(imsalout_index,crit_0);
					if (crit_0>maxsaliency)
					{
						maxsaliency=crit_0;
					}
					//crit_1=votee_tensor.saliency(1);
					//crit_2=votee_tensor.saliency(2);
					dividerv++;
					sums=sums+crit_0;
					//printf("the stick saliency is %f \n", crit_0);
					if(crit_0>thre)
						//if ((crit_0>crit_1)&&(crit_0>crit_2)&&(crit_0>thre))
						//if ((crit_0>crit_1)&&(crit_0>crit_2))
					{
						data Vres;
						Vres.x=cox;
						Vres.y=coy;
						Vres.z=coz;
						queuep.push_back(Vres); 
						C(coz,coy,cox)=1;
					}

				}
				means=sums/dividerv;
				result.clear();

				//if (iter==itert-1)
				//{
				//////write the saliency map for iter-th iteration
				DivideImageFilterType::Pointer divideImageFilter2 = DivideImageFilterType::New();
				divideImageFilter2->SetInput1(imsalout);
				divideImageFilter2->SetInput2(maxsaliency);
				divideImageFilter2->Update();

				typedef itk::ImageFileWriter<ImageTypeFloat> FileWriterTypeFloat;
				FileWriterTypeFloat::Pointer floatwriter2 = FileWriterTypeFloat::New();
				floatwriter2->SetFileName(saliencymap);
				floatwriter2->SetInput(divideImageFilter2->GetOutput());
				try
				{
					floatwriter2->Update();
				}
				catch (itk::ExceptionObject &e)
				{
					std::cout << e << std::endl;
					return 1;
				}
				//}

				//////write the segmentation result for iter-th iteration
				printf("About to print the data for %d-th iteration\n",iter);
				ImageType::Pointer imout2 = ImageType::New();
				imout2->SetRegions(im->GetLargestPossibleRegion());
				imout2->Allocate();
				ImageType::IndexType output_index2;
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
										output_index2[0] = xl+cx-1;
										output_index2[1] = yl+cy-1;
										output_index2[2] = zl+coz-window1;

										//fprintf(fp,"%d %d %d %d",zl+coz+1-window1,yl+cy,xl+cx,M(coz,coy,cox));
										imout2->SetPixel(output_index2,255);
										//#ifdef COMPUTE_L2
										//fprintf(fp," %lf",L(coz,coy,cox));
										//#endif
										//fprintf(fp,"\n");
									}

							}
						}
						char filena_out_iter[1024],iterna[1024];
						strcpy(filena_out_iter,argv[2]);
						sprintf(iterna,"_%dth_iteration_output.tif",iter);
						strcpy(&filena_out_iter[strlen(filena_out_iter)-4],iterna);

						typedef itk::ImageFileWriter<ImageType> FileWriterType;
						FileWriterType::Pointer writer2 = FileWriterType::New();
						writer2->SetFileName(filena_out_iter);
						writer2->SetInput(imout2);
						try
						{
							writer2->Update();
						}
						catch (itk::ExceptionObject &e)
						{
							std::cout << e << std::endl;
							return 1;
						}

				Vm2.clear();
				//free(Cvotee);
				printf("the average stick saliency is %f \n", means);
				printf("%d-th iteration complete\n",iter);
				printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
				time_cal[iter+1]=time(NULL)-t1;

			}
			///////----------------------------------------end of tensor voting------------------------------------------------


			//// Decompose the result for Densification.
			//int rescounter=0;
			//double sums=0;
			//double means;
			//double dividerv=0;
			//double crit_0,crit_1,crit_2;
			//for(int coz = window1;coz<pdepth-window1;coz++)
			//	for(int coy = window; coy<pheight-window;coy+=1)
			//		for(int cox = window; cox<pwidth-window;cox+=1)
			//		{
			//			if (C(coz,coy,cox)==1)
			//			{
			//				rescounter++;
			//				continue;
			//			}
			//			votee_matrix=Vm2[rescounter].voting_matrix;
			//			rescounter++;
			//			rtvl_tensor<3> votee_tensor(votee_matrix);
			//			crit_0=votee_tensor.saliency(0);
			//			crit_1=votee_tensor.saliency(1);
			//			crit_2=votee_tensor.saliency(2);
			//			dividerv++;
			//			sums=sums+crit_0;
			//			printf("the stick saliency is %f \n", crit);
			//			//if(crit_0>thre)
			//			//if ((crit_0>crit_1)&&(crit_0>crit_2)&&(crit_0>thre))
			//			if ((crit_0>crit_1)&&(crit_0>crit_2))
			//			{
			//				data Vres;
			//				Vres.x=cox;
			//				Vres.y=coy;
			//				Vres.z=coz;
			//				queuep.push_back(Vres);
			//				C(coz,coy,cox)=1;
			//			}

			//		}
			//		means=sums/dividerv;
			//		printf("the average stick saliency is %f \n", means);



			///////////////////////////////////////////////////////////////////////////
			//
			//end of tensor voting

			//printf("About to print the data\n");
			//ImageType::Pointer imout = ImageType::New();
			//imout->SetRegions(im->GetLargestPossibleRegion());
			//imout->Allocate();
			//ImageType::IndexType output_index;
			//for(int coz = window1;coz<pdepth-window1;coz++)
			//	for(int coy = window; coy<pheight-window;coy+=1)// skip is used for sampling sometimes. We dont sample in Z direction
			//		for(int cox = window; cox<pwidth-window;cox+=1)
			//		{
			//			if(C(coz,coy,cox)==1)
			//			{
			//				//int minimum_x=MAX(1,cox-window),minimum_y=MAX(1,coy-window),minimum_z=MAX(1,coz-window1);
			//				//int maximum_x=MIN(rwidth/gridsize_x,cox+1-window),maximum_y=MIN(rlength/gridsize_y,coy+1-window),maximum_z=MIN(rdepth,coz+1-window1);

			//				for(int cx = cox+1-window;cx<=cox+1-window;cx++)
			//					for(int cy = coy+1-window;cy<=coy+1-window;cy++)
			//					{
			//						output_index[0] = xl+cx-1;
			//						output_index[1] = yl+cy-1;
			//						output_index[2] = zl+coz-window1;

			//						fprintf(fp,"%d %d %d %d",zl+coz+1-window1,yl+cy,xl+cx,M(coz,coy,cox));
			//						imout->SetPixel(output_index,255);
			//						//#ifdef COMPUTE_L2
			//						fprintf(fp," %lf",L(coz,coy,cox));
			//						//#endif
			//						fprintf(fp,"\n");
			//					}

			//			}
			//		}
			//typedef itk::ImageFileWriter<ImageType> FileWriterType;
			//FileWriterType::Pointer writer = FileWriterType::New();
			//writer->SetFileName(argv[2]);
			//writer->SetInput(imout);
			//try
			//{
			//	writer->Update();
			//}
			//catch (itk::ExceptionObject &e)
			//{
			//	std::cout << e << std::endl;
			//	return 1;
			//}

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
			for (int caltime=0;caltime<itert+1;caltime++)
			{
				printf("%d-th output took %ld seconds\n",caltime,time_cal[caltime]);
			}
			printf("%ld seconds \n",time(NULL)-t1);
			return 0;
}


/////*=========================================================================
////Copyright 2009 Rensselaer Polytechnic Institute
////Licensed under the Apache License, Version 2.0 (the "License");
////you may not use this file except in compliance with the License.
////You may obtain a copy of the License at
////
////http://www.apache.org/licenses/LICENSE-2.0
////
////Unless required by applicable law or agreed to in writing, software
////distributed under the License is distributed on an "AS IS" BASIS,
////WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
////See the License for the specific language governing permissions and
////limitations under the License. 
////=========================================================================*/
//
////#include <stdio.h>
//#include <cstdio>
//#include <vector>
//
//#include <algorithm>
//#include <math.h>
//#include <string.h>
//#include <string>
//#include <sstream>
//#include <assert.h>
//#include <ctime>
//#include "omp.h"
//#include "find_median.cpp"
//
//#include "itkImage.h"
//#include "itkImageFileReader.h"
//#include "itkImageFileWriter.h"
//#include <itkImageRegionIterator.h>
//#include <itkDanielssonDistanceMapImageFilter.h>
//
////lets try tensor 66
//#include <rtvl/rtvl_tensor.hxx>
//#include <rtvl/rtvl_vote.hxx>
//#include <rtvl/rtvl_votee.hxx>
//#include <rtvl/rtvl_voter.hxx>
//#include <rtvl/rtvl_weight_smooth.hxx>
//
////#include <rtvl_tensor.hxx>
////#include <rtvl_vote.hxx>
////#include <rtvl_votee.hxx>
////#include <rtvl_voter.hxx>
////#include <rtvl_weight_smooth.hxx>
//
//
//#include <vnl/vnl_vector_fixed.h>
//#include <vnl/vnl_matrix_fixed.h>
//
//#include <vcl_iostream.h>
////66
//
////for vesselness
//#include "itkHessianToObjectnessMeasureImageFilter.h"
//#include "itkHessianToObjectnessMeasureImageFilter.txx"
//#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
//#include "itkMultiScaleHessianBasedMeasureImageFilter.txx"
//#include "itkCastImageFilter.h"
//#include "itkStatisticsImageFilter.h"
//#include "itkDivideImageFilter.h"
//
//typedef itk::Image<unsigned char,3> ImageType;
//typedef itk::Image<float,3> ImageTypeFloat;
//typedef itk::HessianToObjectnessMeasureImageFilter<float, 3> ObjectnessFilterType;
//typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageTypeFloat, ObjectnessFilterType, ImageTypeFloat> MultiScaleHessianFilterType;
//
//
//
//using namespace std;
//
////#define COMPUTE_L2
//
//#define MAX(a,b) (((a)>(b))?(a):(b))
//#define MIN(a,b) (((a)<(b))?(a):(b))
//
//// The following lines help us access the linear arrays raster and p like a multidimensional array
//#define MAT(a,b,c) raster[(a)*npixels+(b)*w+(c)]
//#define P(a,b,c) p[(a)*pimsize+(b)*pwidth+(c)]
//#define M(a,b,c) matrix[(a)*pimsize+(b)*pwidth+(c)]
//#define C(a,b,c) checked[(a)*pimsize+(b)*pwidth+(c)]
//#define SM(a,b,c) sparmark[(a)*pimsize+(b)*pwidth+(c)]
//#define Mvm1(a,b,c) Vm1[(a)*pimsize+(b)*pwidth+(c)]
//#define Mvm2(a,b,c) Vm2[(a)*pimsize+(b)*pwidth+(c)]
//#define CVOTEE(a,b,c) Cvotee[(a)*pimsize+(b)*pwidth+(c)]
////#define Que0(a,b,c) queuep0.at((a)*pimsize+(b)*pwidth+(c))
////#define MItin(a,b,c)   mark_iterin[(a)*pimsize+(b)*pwidth+(c)]
////#define MItout(a,b,c)  mark_iterout[(a)*pimsize+(b)*pwidth+(c)]
////#ifdef COMPUTE_L2
//#define L(a,b,c) lmatrix[(a)*pimsize+(b)*pwidth+(c)]
////#endif
//
//// window - max window width for X,Y directions
//int window=5,window_fix=5;
//// window1 - max window width for Z direction
//int window1=3;
//// lower threshold
//double lambda_lower_threshold=11,lambda_lower_threshold_fix = 11;
//// higher threshold
//double lambda_higher_threshold=14,lambda_higher_threshold_fix = 14;
////prune value
//double prune=3,prune_adaptiveL=12,prune_adaptiveH=120;
//double epsilonw=1.5,alpha1;
/////threshold for analyse decomposing results
//double thre = 0.8;//0.4//0.5//0.6//0.8 best//1.0//1.8//2.8//2.4//2.2//<1.6 lead to complete
////sparsize for sparse voting
//int sparsize=5,sparsize1=2;///10_2,5_2is the same but much quicker
//const int itert=8;//set iteration time
//#define CFH(a,b) (((unsigned int)a)<<(b*4))
//
//int psize,pimsize,pwidth ,pheight,pdepth ;
//double planes[42][3];
//unsigned char * null,*fore,*back;
//unsigned char *p;
//#define DEBUG_ printf
////compare function for sorting the data
//int compare(const void *a,const void *b)
//{
//	return (int)(*(unsigned char*)a)-(int)(*(unsigned char*)b);
//}
//int compare_double(const void *a,const void *b)
//{
//	return (*(double*)a - *(double*)b);
//}
//struct data{
//	int x,y,z;
//};
//
////struct hypo{
////	vector<double> t;
////	vector<unsigned char> null,back,fore;
////	double amax,bmax,cmax,a,b,c;
////	double L1,likelihood,best_t1,best_t0;
////};
//
//struct vmatrix{
//	vnl_matrix_fixed<double, 3, 3> voting_matrix;
//};
//
//struct return_data{
//	int M;
//	double lvalue;
//};
//void ParseArguments(int argc, char **argv)
//{
//	for (int counter=3; counter<argc; counter++)
//	{
//		istringstream s(argv[counter]);
//		//char name[100];
//		char ch;
//		int d;
//		sscanf(argv[counter],"%c=%d",&ch,&d);
//		switch(ch)
//		{
//		case 'x': case 'y': window = d;
//			DEBUG_("window = %d\n",d);
//			break;
//		case 'z':	window1 = d;
//			DEBUG_("window1 = %d\n",d);
//			break;
//		case 'l':	lambda_lower_threshold =d;
//			DEBUG_("lambda_lower_threshold = %d\n",d);
//			break;
//		case 'h':	lambda_higher_threshold=d;
//			DEBUG_("lambda_higher_threshold = %d\n",d);
//			break;
//		case 'p':	prune=d;
//			DEBUG_("prune = %d\n",d);
//			break;
//		case 'w': epsilonw=d/2.0;
//			DEBUG_("epsilonw = %lf\n",d/2.0);
//			break;
//		default:
//			printf("Parse error in %s\n",argv[counter]);
//			return;
//		}
//	}
//	printf("\n\nSuccessfully parsed %d parameters\n\n\n",argc-1);
//}
//
//vector<data> sparse_mark(int sparsize,int sparsize1,bool * mark_input,bool * mark_output,vector<data> queuep,vector<data> que_output)
//{
//	int poz,poy,pox;
//	printf("Start Neighbor Region Mark for Sparse Voting\n");
//	printf("the input queue size is %d\n",queuep.size());
//	que_output.clear();
//	//vector <data> que_output(queuep);
//	//que_output=queuep;
//	//que_output.assign(queuep.begin(),queuep.end());
//	//printf("que_output size is %d\n",que_output.size());
//	bool * sparmark;
//	sparmark = (bool *) malloc(psize*sizeof(bool));
//	if (sparmark==NULL)
//		printf("Fatal problem: could not malloc enough memory for sparmark\n");
//	memcpy(sparmark,mark_input,psize*sizeof(bool));
//	int control_spar=queuep.size();
//	//printf("que_output size is %d\n",que_output.size());
//
//
//
//	for (int qcounter=0;qcounter<control_spar;qcounter++)
//	{
//		poz=queuep[qcounter].z;
//		poy=queuep[qcounter].y;
//		pox=queuep[qcounter].x;
//		/*if ((poz<window1)||(poy<window)||(pox<window))
//			continue;*/
//		for (int cz=MAX(poz-sparsize1,window1);cz<MIN(poz+sparsize1,pdepth-window1);cz++)
//			for (int cy=MAX(poy-sparsize,window);cy<MIN(poy+sparsize,pheight-window);cy++)
//				for (int cx=MAX(pox-sparsize,window);cx<MIN(pox+sparsize,pwidth-window);cx++)
//				{
//					if (SM(cz,cy,cx)==0)
//					{
//						data po;
//						po.z=cz;
//						po.y=cy;
//						po.x=cx;
//						que_output.push_back(po);
//						SM(cz,cy,cx)=1;
//					}
//
//				}
//
//	}
//	memcpy(mark_output,sparmark,psize*sizeof(bool));
//	free(sparmark);
//	printf("sparse mark complete\n");
//	printf("que_output size is %d\n",que_output.size());
//	return que_output;
//
//
//}
//
//void set_limits(int &zl, int &zu, int &yl, int &yu, int &xl, int&xu,int argc, char**argv)
//{
//	zl = atoi(argv[argc-6]);
//	zu = atoi(argv[argc-5]);
//	yl = atoi(argv[argc-4]);
//	yu = atoi(argv[argc-3]);
//	xl = atoi(argv[argc-2]);
//	xu = atoi(argv[argc-1]);
//
//}
//int main(int argc, char**argv)
//{
////	time_t t1;
////	time_t time_cal[itert+1];
////	t1 = time(NULL);
////	unsigned char * raster;
////	//unsigned char *maxintensity;
////	unsigned int rwidth,rlength,rdepth;
////	// pic file name
////
////	printf("%d\n",argc);
////	//if(argc >2)
////	ParseArguments(argc,argv);
////	//char picfile[1024];// dont segfault me plz!
////	//printf("%s\n",argv[1]);
////	char filenamebuff[1024];//="vessel.pic";
//////	sprintf(filenamebuff,"%s.pic",argv[1]);
////	typedef itk::Image<unsigned char,3> ImageType;
////	typedef itk::ImageFileReader<ImageType> FileReaderType;
////
////	FileReaderType::Pointer reader = FileReaderType::New();
////	reader->SetFileName(argv[1]);
////	try
////	{
////		reader->Update();
////	}
////	catch (itk::ExceptionObject &e)
////	{
////		std::cout << e << std::endl;
////		return 1;
////	}
////
////	ImageType::Pointer im = reader->GetOutput();
////
////	std::cout << "Image Loaded";
////
////	rwidth = im->GetLargestPossibleRegion().GetSize()[0];
////	rlength = im->GetLargestPossibleRegion().GetSize()[1];
////	rdepth = im->GetLargestPossibleRegion().GetSize()[2];
////
////	//itk::ImageRegionIterator<ImageType> imit(im,im->GetRequestedRegion());
////	//#pragma omp parallel for num_threads(14)
////	//for (imit=imit.Begin();imit<=imit.End();++imit)
////	//{
////	//}
////
////
////
////	////for vesselness
////	///////////////////////////////////////
////	//typedef itk::CastImageFilter< ImageType, ImageTypeFloat > CastFilterType;
////	//CastFilterType::Pointer castFilter = CastFilterType::New();
////	//castFilter->SetInput(im);
////	//castFilter->Update();
////	//ImageTypeFloat::Pointer floatout=castFilter->GetOutput();
////	//
////	//typedef itk::StatisticsImageFilter<ImageTypeFloat> StatisticsImageFilterType;
////	//StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New ();
////	//statisticsImageFilter->SetInput(floatout);
////	//statisticsImageFilter->Update();
////	//float maximus=statisticsImageFilter->GetMaximum();
////
////
////
////	//float sigma_min = 1.0f; //0.5f;
////	//float sigma_max = 10.0f; //4.0f;
////	//int sigma_intervals=10;
////	////int sigma_steps = 5;
////
////	//float alpha = 0.5, beta = 0.5, gamma = 0.25; //5.0//0.025;
////	//alpha=alpha*maximus;
////	//beta=beta*maximus;
////	//gamma=gamma*maximus;
////
////	////int obj_dim = objectness_type; //1; //0: Blobness, 1: Vesselness, 2: Plateness
////
////	//MultiScaleHessianFilterType::Pointer multi_scale_Hessian = MultiScaleHessianFilterType::New();
////	//multi_scale_Hessian->SetInput(floatout);
////	//multi_scale_Hessian->SetSigmaMin(sigma_min);
////	//multi_scale_Hessian->SetSigmaMax(sigma_max);
////	//multi_scale_Hessian->SetNumberOfSigmaSteps(sigma_intervals);
////
////	////ObjectnessFilterType::Pointer objectness_filter = ObjectnessFilterType::New();
////	//ObjectnessFilterType::Pointer objectness_filter = multi_scale_Hessian->GetHessianToMeasureFilter();
////	//objectness_filter->SetScaleObjectnessMeasure(false);
////	//objectness_filter->SetBrightObject(true);
////	//objectness_filter->SetAlpha(alpha);
////	//objectness_filter->SetBeta(beta);
////	//objectness_filter->SetGamma(gamma);
////	//objectness_filter->SetObjectDimension(1);/////);
////	////std::cout << obj_measures.alpha << std::endl << obj_measures.beta << std::endl << obj_measures.gamma << std::endl;
////
////	//multi_scale_Hessian->Update();
////
////	//ImageTypeFloat::Pointer floatptr = multi_scale_Hessian->GetOutput();
////
////	//StatisticsImageFilterType::Pointer statisticsImageFilter2 = StatisticsImageFilterType::New ();
////	//statisticsImageFilter2->SetInput(floatptr);
////	//statisticsImageFilter2->Update();
////	//float maximus2=statisticsImageFilter2->GetMaximum();
////	//
////	////itk::Index<3> idx;
////	////idx[0]=55;
////	////idx[1]=68;
////	////idx[2]=21;
////	////float heihei=floatptr->GetPixel(idx);
////	////std::cout<<heihei<<std::endl;
////	typedef itk::DivideImageFilter <ImageTypeFloat, ImageTypeFloat, ImageTypeFloat > DivideImageFilterType;
////	//DivideImageFilterType::Pointer divideImageFilter = DivideImageFilterType::New();
////	//divideImageFilter->SetInput1(floatptr);
////	//divideImageFilter->SetInput2(maximus2);
////	//divideImageFilter->Update();
////
////
////
////	//char vesselnessout[1024];
////	//sprintf(vesselnessout,"%s_vesselness.mhd",argv[1]);
////	//typedef itk::ImageFileWriter<ImageTypeFloat> FileWriterTypeFloat;
////	//FileWriterTypeFloat::Pointer floatwriter = FileWriterTypeFloat::New();
////	//floatwriter->SetFileName(vesselnessout);
////	//floatwriter->SetInput(divideImageFilter->GetOutput());
////	//try
////	//{
////	//	floatwriter->Update();
////	//}
////	//catch (itk::ExceptionObject &e)
////	//{
////	//	std::cout << e << std::endl;
////	//	return 1;
////	//}
////	///////////////////////////////////////////////
////
////	/*
////	FILE * fpi = fopen(filenamebuff,"rb");
////	if(fpi == NULL)
////	{
////		printf("I tried this file '%s'\n",filenamebuff);
////		printf("Couldn't open the file. Please check the filename and try again\n");
////		return 0;
////	}
////
////	unsigned char *temp=new unsigned char[77];
////	fread(temp,sizeof(unsigned char),76,fpi);
////
////	rwidth = CFH(temp[0],0)+CFH(temp[1],2);
////	rlength = CFH(temp[2],0)+CFH(temp[3],2);
////	rdepth = CFH(temp[4],0)+CFH(temp[5],2);
////	*/
////
////	printf("I found rdepth, rlength, rwidth %d %d %d\n",rdepth,rlength,rwidth);
////	//scanf("%*d");
////	int npixels = rwidth*rlength;
////	int w = rwidth;
////	raster = (unsigned char*)malloc(npixels*rdepth*sizeof(unsigned char));
////	if(raster==NULL)
////	{
////		printf("memory problem: couldn't allocate enough memory\n");
////		_exit(0);
////	}
////	memcpy(raster,im->GetBufferPointer(),npixels*rdepth*sizeof(unsigned char));
////
////
////	char buff [1024];// = "E:\\Arun\\My matlab codes\\deconvolved\\deconv";
////	strcpy(buff,filenamebuff);
////
////
////	const int skip  = 1;
////	char filenameoutbuff[1024];
////	sprintf(filenameoutbuff,"%s.npts",argv[1]);
////	FILE * fp = fopen(filenameoutbuff,"w");
////	
////	double *lmatrix;
////	int * matrix;
////	bool *checked;
////	double *l1matrix;
////	double *l2matrix;
////	const int gridsize_x = 1;
////	const int gridsize_y = 1;
////	const int gridx = 0;
////	const int gridy = 0;
////	alpha1 = 0.3;
////
////	vector <data> queuep;
////	//vector <data> queuep0;
////	queuep.reserve(200000);
////	int neighbors[26][3] = {
////		{-1,    -1,    -1},
////		{	-1,     0,    -1},
////		{	-1,     1,    -1},
////		{	0,    -1,    -1},
////		{	0,     0,    -1},
////		{	0,     1,    -1},
////		{	1,    -1,    -1},
////		{	1,     0,    -1},
////		{	1,     1,    -1},
////		{	-1,    -1,     0},
////		{	-1,     0,     0},
////		{	-1,     1,     0},
////		{	0,    -1,     0},
////		{	0,     1,     0},
////		{	1,    -1,     0},
////		{	1,     0,     0},
////		{	1,     1,     0},
////		{	-1,    -1,     1},
////		{	-1,     0,     1},
////		{	-1,     1,     1},
////		{	0,    -1,     1},
////		{	0,     0,     1},
////		{	0,     1,     1},
////		{	1,    -1,     1},
////		{	1,     0,     1},
////		{	1,     1,     1}
////	};
////
////	int zl = 0,zu=rdepth-1,yl=gridy*(rlength/gridsize_y),yu=(gridy+1)*(rlength/gridsize_y)-1,xl=gridx*(rwidth/gridsize_x),xu=(gridx+1)*(rwidth/gridsize_x)-1;
////
////	printf("zl zu yl yu xl xu = %d %d %d %d %d %d\n",zl,zu,yl,yu,xl,xu);
////
////	// We find the various dimensions of the image for our ease of use
////	psize = (zu-zl+1+2*window1)*(yu-yl+1+2*window)*(xu-xl+1+2*window);
////	pimsize = (yu-yl+1+2*window)*(xu-xl+1+2*window);
////	pwidth = (xu-xl+1+2*window);
////	pheight = (yu-yl+1+2*window);
////	pdepth = (zu-zl+1+2*window1);
////	int gamma_size = (2*window+1)*(2*window1+1)*(2*window+1);
////	// allocate memory and copy the data into a new array
////	p = (unsigned char *)malloc(psize*sizeof(unsigned char));
////	matrix = (int *) malloc(psize*sizeof(int));
////	memset(matrix,0,psize*sizeof(int));
////	printf("did I cross that?\n");
////	checked = (bool *) malloc(psize*sizeof(bool));
////	l1matrix = (double *)malloc(gamma_size*sizeof(double));
////	l2matrix = (double *)malloc(gamma_size*sizeof(double));
////	//#ifdef COMPUTE_L2
////	lmatrix = (double *)malloc(psize*sizeof(double));
////	if(l1matrix == NULL || l2matrix == NULL || lmatrix == NULL || matrix == NULL)
////	{
////		printf("mem alloc problem %p %p %p %p\n",l1matrix,l2matrix,lmatrix,matrix);
////
////	}
////	else
////	{
////		printf("No mem alloc problem\n");
////	}
////	//#endif
////	for(int coz = 0; coz <pdepth; coz++)
////		for(int coy = 0 ; coy <pheight; coy++)
////			for(int cox = 0; cox <pwidth; cox++)
////				M(coz,coy,cox)=0;
////
////	printf("Problem here? pwidth = %d pheight = %d pdepth = %d psize = %d \n",pwidth, pheight, pdepth,psize);
////	// replicated
////	int count = 0;
////	for (int coz=window1;coz<pdepth-window1;coz++)
////		for(int coy=window;coy<pheight-window;coy++)
////			for(int cox=window;cox<pwidth-window;cox++)
////			{
////				if(MAT(zl+coz-window1,yl+coy-window,xl+cox-window)>0)
////				{
////					count++;
////				}
////				P(coz,coy,cox)=MAT(zl+coz-window1,yl+coy-window,xl+cox-window);
////			}
////
////			printf("Problem here? pwidth = %d pheight = %d pdepth = %d\n",pwidth, pheight, pdepth);
////			//'replicate' the boundaries
////			for( int counter = 0; counter < window; counter++)
////			{
////				for(int counter1 =0; counter1 < pdepth; counter1++)
////					for(int counter2 = 0; counter2 < pwidth; counter2++)
////					{
////						P(counter1,counter,counter2)=P(counter1,window,counter2);
////						P(counter1,pheight-counter-1,counter2)=P(counter1,pheight-window-1,counter2);
////					}
////					for(int counter1 =0; counter1 < pdepth; counter1++)
////						for(int counter2 = 0; counter2 < pheight; counter2++)
////						{
////							P(counter1,counter2,counter)=P(counter1,counter2,window);
////							P(counter1,counter2,pwidth-counter-1)=P(counter1,counter2,pwidth-window-1);
////						}
////			}
////			for( int counter = 0; counter < window1; counter++)
////			{
////				for(int counter1 =0; counter1 < pheight; counter1++)
////					for(int counter2 = 0; counter2 < pwidth; counter2++)
////					{
////						P(counter,counter1,counter2)=P(window1,counter1,counter2);
////						P(pdepth-counter-1,counter1,counter2)=P(pdepth-window-1,counter1,counter2);
////					}
////			}
////			printf("mat >0 in %d points\n",count);
////			printf("Did I come till epsilon?\n");
////			// continue adding more variables
////			//						double alpha1 = 0.3;
////			//double epsilon = 1e-2;
////
////			//double scaling_matrix[] = {window1, window, window};
////
////			double planes_21[21][3]= {
////        {1, 0, -1},
////				{1, 0, 1},
////				{1, -1, 0},
////				{1, 1, 0},
////				{0, 1, 1},
////				{0, -1, 1},
////				{1, 0, 0},
////				{0, 1, 0},
////				{0, 0, 1},
////				{0.5, 0.5, 1},
////				{0.5, -0.5, 1},
////				{-0.5, -0.5, 1},
////				{-0.5, 0.5, 1},
////				{1, 0.5, 0.5},
////				{1, 0.5, -0.5},
////				{1, -0.5, 0.5},
////				{1, -0.5, -0.5},
////				{0.5, 1, 0.5},
////				{-0.5, 1, 0.5},
////				{0.5, 1, -0.5},
////				{-0.5, 1, -0.5}
////        };
////
////			// Reflect and store in a new array
////			for(int counter = 0; counter < 21; counter++)
////			{
////				planes[counter][0]=planes_21[counter][0];
////				planes[counter][1]=planes_21[counter][1];
////				planes[counter][2]=planes_21[counter][2];
////				planes[counter+21][0]=-planes_21[counter][0];
////				planes[counter+21][1]=-planes_21[counter][1];
////				planes[counter+21][2]=-planes_21[counter][2];
////			}
////			for(int counter =0; counter < 21; counter++)
////			{
////				double sum  = sqrt(planes[counter][0]*planes[counter][0]+planes[counter][1]*planes[counter][1]+planes[counter][2]*planes[counter][2]);
////				planes[counter][0]/=sum;
////				planes[counter][1]/=sum;
////				planes[counter][2]/=sum;
////			}
////			//unsigned int gamma_size = (2*window1+1)*(2*window+1)*(2*window+1);
////			//double t[3]; // stores the lambda values
////
////			/*null = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
////			fore = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
////			back = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));*/
//////			like = (double *)malloc((2*window1+1)*(2*window+1)*(2*window+1)*sizeof(float));
////			//epsilonw=2.5;
////			printf("epsilongw = %lf\n",epsilonw);
////			queuep.clear();
////			double best_t1,best_t0;
////			//double likelihood=0;
////			int NeighborVectorLength=(2*window+1)*(2*window+1)*(2*window1+1);
////			//int StructVectorLength=(pheight-2*window)*(pwidth-2*window)*(pdepth-2*window1);
////			//vector<hypo> hypotest;
////			//hypotest.resize(StructVectorLength);
////			//queuep0.resize(Que0Length);
////			/*vector<unsigned char> null,fore,back;
////			null.resize(NeighborVectorLength);
////			fore.resize(NeighborVectorLength);
////			back.resize(NeighborVectorLength);
////			vector<double> t;
////			t.resize(3);*/
////
////            #pragma omp parallel for num_threads(14)
////			//#pragma omp parallel for 16
////			//#pragma omp parallel for num_threads(10) schedule(dynamic, 1)
////			for(int coz = window1;coz<pdepth-window1;++coz)//coz+=skip
////			{
////				//ImageType::SizeType regionSize;
////				//regionSize[0]=pwidth-2*window;
////				//regionSize[1]=pheight-2*window;
////				//regionSize[2]=1;
////
////				//ImageType::IndexType regionIndex;
////				//regionIndex[0]=window;
////				//regionIndex[1]=window;
////				//regionIndex[2]=coz;
////
////				//ImageType::RegionType region;
////				//region.SetSize(regionSize);
////				//region.SetIndex(regionIndex);
////
////				//itk::ImageRegionIterator<ImageType> imit(im,region);
////				//#pragma omp parallel for num_threads(10) schedule(dynamic, 1)
////				//for (imit = imit.Begin(); !imit.IsAtEnd(); ++imit)
////				for(int coy = window; coy<pheight-window;++coy)//coy+=skip skip is used for sampling. We dont sample in Z direction
////					for(int cox = window; cox<pwidth-window;++cox)//cox+=skip
////					{
////						//hypo hy;
////						//double t[3];
////						vector<double> t;
////						t.resize(3);
////						/*null = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
////						fore = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));
////						back = (unsigned char*)malloc((2*window1+1)*(2*window+1)*(2*window+1));*/
////						vector<unsigned char> null,fore,back;
////						null.resize(NeighborVectorLength);
////						fore.resize(NeighborVectorLength);
////						back.resize(NeighborVectorLength);
////						/*unsigned char *null = 0;
////						null = new unsigned char[];*/
////						/*unsigned char *null=new unsigned char[847];
////						unsigned char *back=new unsigned char[847];
////						unsigned char *fore=new unsigned char[847];*/
////						double best_t1,best_t0;
////						best_t1 = best_t0 = -1;
////						//double t[3];
////
////						//the main for loop
////						double L1 = -5000000,likelihood=0;
////						//double L;
////						double amax,bmax,cmax;
////						//int pcount =0;
////						int pccn =0;
////
////
////						//ImageType::SizeType regionSize2;
////						//regionSize2[0]=2*window+1;
////						//regionSize2[1]=2*window+1;
////						//regionSize2[2]=2*window1+1;
////
////						//ImageType::IndexType regionIndex2;
////						//regionIndex2[0] = imit.GetIndex()[0]-window;
////						//regionIndex2[1] = imit.GetIndex()[1]-window;
////						//regionIndex2[2] = imit.GetIndex()[2]-window1;
////
////						//ImageType::RegionType region2;
////						//region.SetSize(regionSize2);
////						//region.SetIndex(regionIndex2);
////
////						//itk::ImageRegionIterator<ImageType> imit2(im,region2);
////
////						//for (imit2=imit2.Begin();!imit2.IsAtEnd();++imit2)
////						for(int wz=-window1;wz<=window1;wz++)
////							for(int wy=-window;wy<=window;wy++)
////								for(int wx=-window;wx<=window;wx++)
////								{
////									/*null.at(pccn)=P(coz+wz,coy+wy,cox+wx);
////									pccn++;*/
////									null[pccn++]=P(coz+wz,coy+wy,cox+wx);
////									//null[pccn++]=imit2.Get();
////								}
////
////						//defined in stdlib to sort sequential arrays
////						//	qsort(null,pccn,sizeof(unsigned char),compare);
////						// find the median of the null hypothesis
////						/*double sum_sort =0;
////						int sum_sort_num =0;
////						for(int counter_s=alpha1*pccn; counter_s<pccn-alpha1*pccn;counter_s++,sum_sort_num++)
////						sum_sort+=null[counter_s];
////						t[2]=sum_sort/sum_sort_num;
////						*/
////						//if(pccn%2==0)
////						//	t[2]=(null[pccn/2]+null[pccn/2-1])/2.0;
////						//else
////						//	t[2]=null[pccn/2];
////							
////						/*int lengthnull=null.size();
////						unsigned char *nullarray;
////						nullarray=new unsigned char[lengthnull];
////						for (int ii=0;ii<lengthnull;++ii)
////						{
////							nullarray[ii]=null.at(ii);
////						}*/
////								t[2] = median(pccn,null);
////						// a small threshold fixed to accelerate the program. Anything with a median of less than 2 is ignored for DEBUG
////
////						if(t[2]<prune)
////							continue;
////						//if(t[2]<prune_adaptiveL)//7
////						//{
////						//	lambda_lower_threshold=lambda_lower_threshold_fix-5;
////						//	lambda_higher_threshold=lambda_higher_threshold_fix-5;
////
////						//}
////						//if(prune_adaptiveH<t[2])//10
////						//{
////						//	lambda_lower_threshold=lambda_lower_threshold_fix+2;
////						//	lambda_higher_threshold=lambda_higher_threshold_fix+2;
////						//	window=window_fix-2;
////						//	//int pccn =0;
////						//	//for(int wz=-window1;wz<=window1;wz++)
////						//	//	for(int wy=-window;wy<=window;wy++)
////						//	//		for(int wx=-window;wx<=window;wx++)
////						//	//		{
////						//	//			null[pccn++]=P(coz+wz,coy+wy,cox+wx);
////						//	//		}
////
////
////						//}
////
////						double a,b,c;
////						//printf("I did come here\n");
////						// now for the alternate hypothesis
////						for(int co = 0; co < 21;co++)
////						{
////
////							a = planes[co][0];
////							b = planes[co][1];
////							c = planes[co][2];
////							int pccf =0;
////							int pccb =0;
////							
////						//ImageType::SizeType regionSize2;
////						//regionSize2[0]=2*window+1;
////						//regionSize2[1]=2*window+1;
////						//regionSize2[2]=2*window1+1;
////
////						//ImageType::IndexType regionIndex2;
////						//regionIndex2[0]=imit.GetIndex()[0]-window;
////						//regionIndex2[1]=imit.GetIndex()[1]-window;
////						//regionIndex2[2]=imit.GetIndex()[2]-window1;
////
////						//ImageType::RegionType region2;
////						//region.SetSize(regionSize2);
////						//region.SetIndex(regionIndex2);
////
////						//itk::ImageRegionIterator<ImageType> imit2(im,region2);
////
////						//for (imit2=imit2.Begin();!imit2.IsAtEnd();++imit2)
////							for(int wz =-window1;wz<=window1;wz++)
////								for(int wy=-window;wy<=window;wy++)
////									for(int wx=-window;wx<=window;wx++)
////									{
////										//ImageType::IndexType neighborIndex;
////										//ImageType::IndexType neighborIndex2;
////										//neighborIndex=imit.GetIndex();
////										//neighborIndex2=imit2.GetIndex();
////										//int wz,wy,wx;
////										//wz=neighborIndex2[2]-neighborIndex[2];
////										//wy=neighborIndex2[1]-neighborIndex[1];
////										//wx=neighborIndex2[0]-neighborIndex[0];
////										//FIXME
////										if((wz*a+wy*b+wx*c)<= epsilonw  && (wz*a+wy*b+wx*c)>=-epsilonw)
////										{
////											/*fore.at(pccf)=P(coz+wz,coy+wy,cox+wx);
////											pccn++;*/
////											fore[pccf++]=P(coz+wz,coy+wy,cox+wx);
////											//fore[pccf++]=imit2.Get();
////										}
////										else
////										{
////											/*back.at(pccb)=P(coz+wz,coy+wy,cox+wx);
////											pccb++;*/
////											back[pccb++]=P(coz+wz,coy+wy,cox+wx);
////											//back[pccb++]=imit2.Get();
////										}
////									}
////									
////
////									//this is an inefficient way to find the median. 
////									//I sort the arrays and find the middle element
////									//qsort(fore,pccf,sizeof(unsigned char),compare);
////									//qsort(back,pccb,sizeof(unsigned char),compare);
////									//find the median of back
////									//if(pccb%2==0)
////									//	t[0]=(back[pccb/2]+back[pccb/2-1])/2.0;
////									//else
////									//	t[0]=back[pccb/2];
////									////find the median of fore
////									//if(pccf%2==0)
////									//	t[1]=(fore[pccf/2]+fore[pccf/2-1])/2.0;
////									//else
////									//	t[1]=fore[pccf/2];
////
////									/*int lengthback=back.size();
////									unsigned char *backarray;
////									backarray=new unsigned char[lengthback];
////									for (int ii=0;ii<lengthback;++ii)
////									{
////										backarray[ii]=back.at(ii);
////									}
////									int lengthfore=fore.size();
////									unsigned char *forearray;
////									forearray=new unsigned char[lengthfore];
////									for (int ii=0;ii<lengthfore;++ii)
////									{
////										forearray[ii]=fore.at(ii);
////									}*/
////									t[0] = median(pccb,back);
////									t[1] = median(pccf,fore);
////								
////									likelihood=0;
////									if(t[1]-t[0]>=lambda_lower_threshold)
////									{
////										//L=t[1]-t[0];
////
////										if(t[2]==0)
////										{
////											likelihood = 1000000;//set some high value because likelihood is \inf
////										}
////										else
////										{
////											//new 
////											likelihood += (t[1]+t[0]-2*t[2])/t[2]+log(t[2]/t[1]);
////											if(t[0]!=0)
////											{
////												likelihood +=log(t[2]/t[0]);
////											}
////											if(likelihood<0)
////											{
////												likelihood = 1000000;
////											}
////										}
////									}
////									else
////										likelihood=-5000000;
////									if(likelihood>L1)
////									{
////										amax = a;
////										bmax = b;
////										cmax = c;
////										best_t1 = t[1];
////										best_t0 = t[0];
////										L1 = likelihood;
////										//if(L1 > 0)
////										//	break;
////									}
////						}
////						//int cox,coy;
////						//ImageType::IndexType regionIndex3;
////						//regionIndex3=imit.GetIndex();
////						//cox=regionIndex3[0];
////						//coy=regionIndex3[1];
////						C(coz,coy,cox)=false;
////						if(L1>0) // tau = 0
////						{
////							//printf("I came here\n");
////							M(coz,coy,cox) = int(best_t1-best_t0+0.5);
////							if(best_t1-best_t0 >= lambda_higher_threshold)
////							{
////								C(coz,coy,cox)=true;
////								//data d ;
////								//d.x = cox;
////								//d.y = coy;
////								//d.z = coz;
////								//#pragma omp critical
////								//{
////								//	//queuep.push_back(d);
////								//	Que0(coz,coy,cox)=d;
////								//}
////							}
////							L(coz,coy,cox)=L1;
////							//compute L2 for the point
////							//double L2 = 0;
////
////						}
////						else
////						{
////							M(coz,coy,cox)=-1;
////
////
////						}
////						//lambda_lower_threshold=lambda_lower_threshold_fix;
////						//lambda_higher_threshold=lambda_higher_threshold_fix;
////						//window=window_fix;
////						/*delete null;
////						delete fore;
////						delete back;*/
////						null.clear();
////						fore.clear();
////						back.clear();
////						t.clear();
////						/*free(null);
////						free(back);
////						free(fore);*/
////					}
////					printf("%s: %d coz\n",filenamebuff,coz);
////					printf("%s: %ld seconds have elapsed since starting\n",filenamebuff, time(NULL)-t1);
////					fflush(stdout);//need revise for omp
////			}
////			// now start with seeds and expand at all places again
////			// Things TODO :
////			// 
////			//free(p);
////			for(int coz = window1;coz<pdepth-window1;++coz)//coz+=skip
////				for(int coy = window; coy<pheight-window;++coy)//coy+=skip skip is used for sampling. We dont sample in Z direction
////					for(int cox = window; cox<pwidth-window;++cox)
////					{
////						if (C(coz,coy,cox)==true)
////						{
////							data d ;
////						    d.x = cox;
////							d.y = coy;
////							d.z = coz;
////							queuep.push_back(d);
////						}
////					}
////			int pc = -1;
////			
////			lambda_lower_threshold=lambda_lower_threshold_fix;
////			while(1)
////			{
////				pc++;
////				if((unsigned int)pc == queuep.size())
////					break;
////				data seedp = queuep[pc];
////				assert(M(seedp.z,seedp.y,seedp.x)>0);
////				int l1 = M(seedp.z,seedp.y,seedp.x);
////
////				//double l2 = L(seedp.z,seedp.y,seedp.x);
////
////				for (int co = 0; co<26; co++) // 3x3x3 neighbors
////				{
////					int x1 = seedp.x+neighbors[co][0];
////					int y1 = seedp.y+neighbors[co][1];
////					int z1 = seedp.z+neighbors[co][2];
////					if(!(x1>=window && y1>=window && z1>=window1 && x1<pwidth-window && y1<pheight-window && z1 <pdepth-window1))
////					{
////						continue;
////					}
////					if(C(z1,y1,x1)==1)
////						continue;
////					
////					if(M(z1,y1,x1)>=MAX(l1-8,lambda_lower_threshold))
////					{
////						data D;
////						D.x = x1;
////						D.y = y1;
////						D.z = z1;
////						queuep.push_back(D);
////						C(z1,y1,x1)=1;
////					}
////				}
////			}
////			printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
////
////			printf("About to print the data after initial segmentation\n");
////			ImageType::Pointer imout = ImageType::New();
////			imout->SetRegions(im->GetLargestPossibleRegion());
////			imout->Allocate();
////			ImageType::IndexType output_index;
////			for(int coz = window1;coz<pdepth-window1;coz++)
////				for(int coy = window; coy<pheight-window;coy+=1)// skip is used for sampling sometimes. We dont sample in Z direction
////					for(int cox = window; cox<pwidth-window;cox+=1)
////					{
////						if(C(coz,coy,cox)==1)
////						{
////							//int minimum_x=MAX(1,cox-window),minimum_y=MAX(1,coy-window),minimum_z=MAX(1,coz-window1);
////							//int maximum_x=MIN(rwidth/gridsize_x,cox+1-window),maximum_y=MIN(rlength/gridsize_y,coy+1-window),maximum_z=MIN(rdepth,coz+1-window1);
////
////							for(int cx = cox+1-window;cx<=cox+1-window;cx++)
////								for(int cy = coy+1-window;cy<=coy+1-window;cy++)
////								{
////									output_index[0] = xl+cx-1;
////									output_index[1] = yl+cy-1;
////									output_index[2] = zl+coz-window1;
////
////									fprintf(fp,"%d %d %d %d",zl+coz+1-window1,yl+cy,xl+cx,M(coz,coy,cox));
////									imout->SetPixel(output_index,255);
////									//#ifdef COMPUTE_L2
////									fprintf(fp," %lf",L(coz,coy,cox));
////									//#endif
////									fprintf(fp,"\n");
////								}
////
////						}
////					}
////			typedef itk::ImageFileWriter<ImageType> FileWriterType;
////			FileWriterType::Pointer writer = FileWriterType::New();
////			writer->SetFileName(argv[2]);
////			writer->SetInput(imout);
////			try
////			{
////				writer->Update();
////			}
////			catch (itk::ExceptionObject &e)
////			{
////				std::cout << e << std::endl;
////				return 1;
////			}
////			printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
////			time_cal[0]=time(NULL)-t1;
////			free(matrix);
////			free(lmatrix);
////			free(l1matrix);
////			free(l2matrix);
////			free(raster);
////			fclose(fp);
//
//			//----------------------function for reading segmented file and use it for tensor voting-------------
//			//-------cannot be with the above part simultaneously
//			time_t t1;
//			time_t time_cal[itert+1];
//			t1 = time(NULL);
//			unsigned char * raster;
//			//unsigned char *maxintensity;
//			unsigned int rwidth,rlength,rdepth;
//			// pic file name
//
//			printf("%d\n",argc);
//			//if(argc >2)
//			ParseArguments(argc,argv);
//			//char picfile[1024];// dont segfault me plz!
//			//printf("%s\n",argv[1]);
//			char filenamebuff[1024];//="vessel.pic";
//			//	sprintf(filenamebuff,"%s.pic",argv[1]);
//			typedef itk::Image<unsigned char,3> ImageType;
//			typedef itk::ImageFileReader<ImageType> FileReaderType;
//
//			FileReaderType::Pointer reader = FileReaderType::New();
//			reader->SetFileName(argv[1]);
//			try
//			{
//				reader->Update();
//			}
//			catch (itk::ExceptionObject &e)
//			{
//				std::cout << e << std::endl;
//				return 1;
//			}
//
//			ImageType::Pointer im = reader->GetOutput();
//
//			std::cout << "Image Loaded"<<endl;
//
//			rwidth = im->GetLargestPossibleRegion().GetSize()[0];
//			rlength = im->GetLargestPossibleRegion().GetSize()[1];
//			rdepth = im->GetLargestPossibleRegion().GetSize()[2];
//
//			/*ImageType::Pointer imout = ImageType::New();
//			imout->SetRegions(im->GetLargestPossibleRegion());
//			imout->Allocate();*/
//
//			typedef itk::DivideImageFilter <ImageTypeFloat, ImageTypeFloat, ImageTypeFloat > DivideImageFilterType;
//
//			vector <data> queuep;
//			queuep.clear();
//			const int gridsize_x = 1;
//			const int gridsize_y = 1;
//			const int gridx = 0;
//			const int gridy = 0;
//			alpha1 = 0.3;
//			int zl = 0,zu=rdepth-1,yl=gridy*(rlength/gridsize_y),yu=(gridy+1)*(rlength/gridsize_y)-1,xl=gridx*(rwidth/gridsize_x),xu=(gridx+1)*(rwidth/gridsize_x)-1;
//
//			psize = (zu-zl+1+2*window1)*(yu-yl+1+2*window)*(xu-xl+1+2*window);
//			pimsize = (yu-yl+1+2*window)*(xu-xl+1+2*window);
//			pwidth = (xu-xl+1+2*window);
//			pheight = (yu-yl+1+2*window);
//			pdepth = (zu-zl+1+2*window1);
//			bool *checked;
//			checked = (bool *) malloc(psize*sizeof(bool));
//			ImageType::IndexType pixelIndex;
//			int pixelval;
//			data ddd;
//			itk::ImageRegionIterator<ImageType> imit(im,im->GetRequestedRegion());
//			for (imit=imit.Begin();imit<=imit.End();++imit)
//			{
//				/*if (imit==imit.Begin())
//				{
//					cout<<"66"<<endl;
//				}*/
//				pixelval=imit.Get();
//				if (pixelval!=0)
//				{
//					
//					ddd.x=imit.GetIndex()[0]+window;
//					ddd.y=imit.GetIndex()[1]+window;
//					ddd.z=imit.GetIndex()[2]+window1;
//					queuep.push_back(ddd);
//					C(ddd.z,ddd.y,ddd.x)=1;
//
//				}
//
//
//			}
//
//
//			//---------------------------------------------------------------------------------------------------
//
//			//-------------------------------------now for tensor voting----------------------------------------------------
//			//
//			///////////////////////////////////////////////////////////////////////
//			
//			for (int iter=0; iter<itert; iter++)
//			{
//				printf("%d-th iteration start...\n",iter);
//				//define ball tensor for initialization
//				vnl_matrix_fixed<double, 3, 3> voter_matrix;
//				//vnl_matrix_fixed<double, 3, 3> votee_matrix(0.0);
//				voter_matrix(0,0) = 1;
//				voter_matrix(0,1) = 0;
//				voter_matrix(0,2) = 0;
//				voter_matrix(1,0) = 0;
//				voter_matrix(1,1) = 1;
//				voter_matrix(1,2) = 0;
//				voter_matrix(2,0) = 0;
//				voter_matrix(2,1) = 0;
//				voter_matrix(2,2) = 1;
//				//vnl_vector_fixed<double,3> voter_location(0.0);
//				//vmatrix temp;
//				int control;
//				control=queuep.size();
//				//int control2;
//				//int StructVectorLength=(pheight-2*window)*(pwidth-2*window)*(pdepth-2*window1)+1;
//				//control=(pwidth-2*window)*(pheight-2*window)*(pdepth-2*window1)-1;
//				vector<vmatrix> Vm1;//store the votee_matrix result for token refinement,corresponding to each point in queuep
//				vector<vmatrix> Vm2;
//				//vector<vmatrix> Vm3;
//				//Vm.reserve(2000000);
//				Vm1.clear();
//				
//				/*std::vector< std::vector < vmatrix > > test;
//				test.resize(control1);*/
//
//				printf("Initialization for tensor voting complete\n");
//
//
//
//				//Token Refinement
//				printf("Start Token Refinement...\n");
//				//control1=queuep.size();
//				Vm1.resize(psize);
//				printf("the queuep size is %d\n",control);
//				//#pragma omp parallel for 16
//				#pragma omp parallel for num_threads(12)
//				for(int coz = window1;coz<pdepth-window1;++coz)//coz+=skip
//					for(int coy = window; coy<pheight-window;++coy)//coy+=skip skip is used for sampling. We dont sample in Z direction
//						for(int cox = window; cox<pwidth-window;++cox)
//				
//						{
//					if (C(coz,coy,cox)!=1)
//						continue;
//					vnl_matrix_fixed<double, 3, 3> votee_matrixaccu1(0.0);
//					vnl_vector_fixed<double, 3> votee_location;
//					votee_location[0]=(double)cox;
//					votee_location[1]=(double)coy;
//					votee_location[2]=(double)coz;
//					for (int cz=MAX(coz-sparsize1,window1);cz<MIN(coz+sparsize1,pdepth-window1);++cz)
//						for (int cy=MAX(coy-sparsize,window);cy<MIN(coy+sparsize,pheight-window);++cy)
//							for (int cx=MAX(cox-sparsize,window);cx<MIN(cox+sparsize,pwidth-window);++cx)//maybe should use different sparsize for token refinement
//							{
//								if (C(cz,cy,cx)!=1)
//									continue;
//								//data voterp = queuep[voterc];
//								vnl_vector_fixed<double,3> voter_location(0.0);
//								voter_location[0]=(double)cx;
//								voter_location[1]=(double)cy;
//								voter_location[2]=(double)cz;
//								vnl_matrix_fixed<double, 3, 3> votee_matrix(0.0);
//								// Use "rtvl_tensor" to decompose the matrix.		
//								rtvl_tensor<3> voter_tensor(voter_matrix);
//
//								// Use "rtvl_voter" to encapsulate a token (location + input tensor).
//								rtvl_voter<3> voter(voter_location, voter_tensor);
//
//								// Use "rtvl_votee" to encapsulate a site (location + output tensor).
//								rtvl_votee<3> votee(votee_location, votee_matrix);
//								//vcl_cout << vcl_endl;
//
//								// Choose a weight profile, initialized with spatial scale.
//								rtvl_weight_smooth<3> tvw(1.0);
//
//								// Compute one vote.
//								rtvl_vote(voter, votee, tvw);
//
//								votee_matrixaccu1+=votee_matrix;
//								//printf("%d-th voter for Token Refinement\n",voterc);
//
//
//							}
//							vmatrix temp;
//							temp.voting_matrix=votee_matrixaccu1;
//							Mvm1(coz,coy,cox)=temp;
//							//Vm1.at(voteec)=temp;
//							//Vm1.push_back(temp);
//							//test.at(voteec).push_back(temp);
//							//printf("%d-th votee complete for Token Refinement\n",voteec);
//				}
//				printf("Token Refinement complete\n");
//				printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
//
//				///Sparse Voting Mark
//
//				//int sparsize=10,sparsize1=2;
//				vector<data> que_output;
//				bool * mark_sparout;
//				mark_sparout = (bool *) malloc(psize*sizeof(bool));
//				if (mark_sparout==NULL)
//					printf("Fatal problem: could not malloc enough memory for mark_sparout\n");
//				que_output=sparse_mark(sparsize,sparsize1,checked,mark_sparout,queuep,que_output);
//				//printf("que_output size is %d\n",que_output.size());
//				printf("Sparse Voting region marked after Token Refinement\n");
//				printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
//
//				///Sparse Voting
//				//int itert=1;
//				//control2=Vm1.size();
//				vector <data> result(que_output);
//				//result=que_output;
//				//vector <data> queue_voter(queuep);
//				//vector <data> queue_votee(que_output);
//				que_output.clear();
//
//				bool * Cvotee;
//				Cvotee = (bool *) malloc(psize*sizeof(bool));
//				if (Cvotee==NULL)
//				printf("Fatal problem: could not malloc enough memory for Cvotee\n");
//				memcpy(Cvotee,mark_sparout,psize*sizeof(bool));
//				free(mark_sparout);
//
//				//bool * mark_iterout;
//				//mark_iterout = (bool *) malloc(psize*sizeof(bool));
//				//if (mark_iterout==NULL)
//				//printf("Fatal problem: could not malloc enough memory for mark_iterout\n");
//
//				printf("Start Sparse Voting...\n");
//				//int voteenum;
//				//vector<data> queueacc;
//				//queueacc.clear();
//				//queueacc=queuep;
//
//				//for (int iter=0; iter<itert; iter++)
//				//{
//				//vector <data> queue_votee(que_output);
//				//printf("%d-th iteration begins...\n",iter);
//				//queue_votee=que_output;
//				//voteenum=queue_votee.size();
//				/*int lenGth;
//				lenGth=queue_votee.size();
//				for (int it=0;it<lenGth;it++)
//				{
//				data DAta=queue_votee[it];
//				queueacc.push_back(DAta);
//
//				}*/
//				Vm2.clear();
//				///////////main loop for votee, one voter to one votee at a time						
//				//int vcounter=0;//store matrix for votee
//				Vm2.resize(psize);
//				//#pragma omp parallel for 16
//				#pragma omp parallel for num_threads(12)
//				for(int coz = window1;coz<pdepth-window1;++coz)//coz+=skip
//					for(int coy = window; coy<pheight-window;++coy)//coy+=skip skip is used for sampling. We dont sample in Z direction
//						for(int cox = window; cox<pwidth-window;++cox)
//				{
//					//printf("did i come here 66_1\n");
//					if (CVOTEE(coz,coy,cox)!=1)
//						continue;
//					vnl_matrix_fixed<double, 3, 3> votee_matrixaccu2(0.0);
//					//vote to one point each time
//					vnl_vector_fixed<double, 3> votee_location;
//					//int tempx,tempy,tempz;
//					//tempx=queue_votee[ivotee].x;
//					//tempy=queue_votee[ivotee].y;
//					//tempz=queue_votee[ivotee].z;
//					votee_location[0] = (double)cox;
//					votee_location[1] = (double)coy;
//					votee_location[2] = (double)coz;
//					//printf("the value of control2 is %d\n",control2);
//
//					/////loop for voter
//					for (int cz=MAX(coz-sparsize1,window1);cz<MIN(coz+sparsize1,pdepth-window1);++cz)
//						for (int cy=MAX(coy-sparsize,window);cy<MIN(coy+sparsize,pheight-window);++cy)
//							for (int cx=MAX(cox-sparsize,window);cx<MIN(cox+sparsize,pwidth-window);++cx)
//					{
//						//printf("did i come here 66_2\n");
//						if (C(cz,cy,cx)!=1)
//							continue;
//						vnl_vector_fixed<double,3> voter_location(0.0);
//						/*int tx,ty,tz;
//						tx=queue_voter[voterc].x;
//						ty=queue_voter[voterc].y;
//						tz=queue_voter[voterc].z;*/
//						voter_location[0]=(double)cx;
//						voter_location[1]=(double)cy;
//						voter_location[2]=(double)cz;
//						/*if (abs(tempx-tx)>sparsize)
//							continue;
//						if (abs(tempy-ty)>sparsize)
//							continue;
//						if (abs(tempz-tz)>sparsize1)
//							continue;*/
//
//						vnl_matrix_fixed<double, 3, 3> votee_matrix(0.0);
//
//
//						//voter_matrix = Vm1[voterc].voting_matrix;
//						voter_matrix = Mvm1(cz,cy,cx).voting_matrix;
//
//						// Use "rtvl_tensor" to decompose the matrix.		
//						rtvl_tensor<3> voter_tensor(voter_matrix);
//
//						// Use "rtvl_voter" to encapsulate a token (location + input tensor).
//						rtvl_voter<3> voter(voter_location, voter_tensor);
//
//						// Use "rtvl_votee" to encapsulate a site (location + output tensor).
//						rtvl_votee<3> votee(votee_location, votee_matrix);
//						//vcl_cout << vcl_endl;
//
//						// Choose a weight profile, initialized with spatial scale.
//						rtvl_weight_smooth<3> tvw(1.0);
//
//						// Compute one vote. 
//						rtvl_vote(voter, votee, tvw);
//
//						votee_matrixaccu2+=votee_matrix;
//						//printf("%d-th voter for tensor",voterc);
//						//printf("did i come here 66_4\n");
//
//					}
//					//printf("did i come here 66_3\n");
//					vmatrix temp;
//					temp.voting_matrix=votee_matrixaccu2;
//					Mvm2(coz,coy,cox)=temp;
//					//Vm2.at(ivotee)=temp;
//					//Vm2.push_back(temp);
//					//vcounter++;
//					//printf("%d-th votee complete for Sparse Voting\n",vcounter);
//
//				}
//				free(Cvotee);
//				//result=queue_votee;
//				//printf("did i come here 66_4\n");
//				//printf("%d-th iteration for tensor",iter);
//				//if (iter==itert-1)
//				//{
//				//	//result=queueacc;
//				//	result=queue_votee;
//				//}
//				//control2=Vm2.size();
//				//vector <data> queue_voter(queue_votee);
//				//queue_voter=queue_votee;
//				//Vm1.clear();
//				//Vm1=Vm2;
//				//call sparse mark function to expand the region
//				//que_output=sparse_mark(sparsize,sparsize1,mark_iterin,mark_iterout,queue_votee,que_output);
//				//memcpy(mark_iterin,mark_iterout,psize*sizeof(bool));				
//				//queue_votee.clear();
//
//				//if (iter<itert-1)
//				// vector <data> queue_votee(que_output);
//				//vector <data> result(queue_votee);
//				//printf("%d-th iteration complete\n",iter);
//				//printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
//
//				//}
//				//free(mark_iterin);
//				//free(mark_iterout);
//				Vm1.clear();
//				printf("Sparse Voting complete\n");
//				printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
//				//Vm2.clear();
//
//				// Decompose the result for sparse voting.
//				printf("Start Decomposing the result for sparse voting\n");
//
//				char saliencymap[1024];
//				sprintf(saliencymap,"%s_saliencymap_%d_iteration.mhd",argv[2],iter);
//				ImageTypeFloat::Pointer imsalout = ImageTypeFloat::New();
//				imsalout->SetRegions(im->GetLargestPossibleRegion());
//				imsalout->Allocate();
//				ImageTypeFloat::IndexType imsalout_index;
//
//				double maxsaliency=0;
//				double sums=0;
//				double means;
//				double dividerv=0;
//				double crit_0,crit_1,crit_2;
//				int voteesize=result.size();
//				int coz,coy,cox;
//				vnl_matrix_fixed<double, 3, 3> votee_matrix(0.0);
//				for (int voteecount=0;voteecount<voteesize;voteecount++)
//				{
//					coz=result[voteecount].z;
//					coy=result[voteecount].y;
//					cox=result[voteecount].x;
//					if (C(coz,coy,cox)==1)
//					{
//						continue;
//					}
//					imsalout_index[0]=cox-window;
//					imsalout_index[1]=coy-window;
//					imsalout_index[2]=coz-window1;
//					//votee_matrix=Vm2[voteecount].voting_matrix;
//					votee_matrix=Mvm2(coz,coy,cox).voting_matrix;
//					rtvl_tensor<3> votee_tensor(votee_matrix);
//					crit_0=votee_tensor.saliency(0);
//					imsalout->SetPixel(imsalout_index,crit_0);
//					if (crit_0>maxsaliency)
//					{
//						maxsaliency=crit_0;
//					}
//					//crit_1=votee_tensor.saliency(1);
//					//crit_2=votee_tensor.saliency(2);
//					dividerv++;
//					sums=sums+crit_0;
//					//printf("the stick saliency is %f \n", crit_0);
//					if(crit_0>thre)
//						//if ((crit_0>crit_1)&&(crit_0>crit_2)&&(crit_0>thre))
//						//if ((crit_0>crit_1)&&(crit_0>crit_2))
//					{
//						data Vres;
//						Vres.x=cox;
//						Vres.y=coy;
//						Vres.z=coz;
//						queuep.push_back(Vres); 
//						C(coz,coy,cox)=1;
//					}
//
//				}
//				means=sums/dividerv;
//				result.clear();
//
//				//if (iter==itert-1)
//				//{
//				//////write the saliency map for iter-th iteration
//				DivideImageFilterType::Pointer divideImageFilter2 = DivideImageFilterType::New();
//				divideImageFilter2->SetInput1(imsalout);
//				divideImageFilter2->SetInput2(maxsaliency);
//				divideImageFilter2->Update();
//
//				typedef itk::ImageFileWriter<ImageTypeFloat> FileWriterTypeFloat;
//				FileWriterTypeFloat::Pointer floatwriter2 = FileWriterTypeFloat::New();
//				floatwriter2->SetFileName(saliencymap);
//				floatwriter2->SetInput(divideImageFilter2->GetOutput());
//				try
//				{
//					floatwriter2->Update();
//				}
//				catch (itk::ExceptionObject &e)
//				{
//					std::cout << e << std::endl;
//					return 1;
//				}
//				//}
//
//				//////write the segmentation result for iter-th iteration
//				printf("About to print the data for %d-th iteration\n",iter);
//				ImageType::Pointer imout2 = ImageType::New();
//				imout2->SetRegions(im->GetLargestPossibleRegion());
//				imout2->Allocate();
//				ImageType::IndexType output_index2;
//				for(int coz = window1;coz<pdepth-window1;coz++)
//					for(int coy = window; coy<pheight-window;coy+=1)// skip is used for sampling sometimes. We dont sample in Z direction
//						for(int cox = window; cox<pwidth-window;cox+=1)
//						{
//							if(C(coz,coy,cox)==1)
//							{
//								//int minimum_x=MAX(1,cox-window),minimum_y=MAX(1,coy-window),minimum_z=MAX(1,coz-window1);
//								//int maximum_x=MIN(rwidth/gridsize_x,cox+1-window),maximum_y=MIN(rlength/gridsize_y,coy+1-window),maximum_z=MIN(rdepth,coz+1-window1);
//
//								for(int cx = cox+1-window;cx<=cox+1-window;cx++)
//									for(int cy = coy+1-window;cy<=coy+1-window;cy++)
//									{
//										output_index2[0] = xl+cx-1;
//										output_index2[1] = yl+cy-1;
//										output_index2[2] = zl+coz-window1;
//
//										//fprintf(fp,"%d %d %d %d",zl+coz+1-window1,yl+cy,xl+cx,M(coz,coy,cox));
//										imout2->SetPixel(output_index2,255);
//										//#ifdef COMPUTE_L2
//										//fprintf(fp," %lf",L(coz,coy,cox));
//										//#endif
//										//fprintf(fp,"\n");
//									}
//
//							}
//						}
//						char filena_out_iter[1024],iterna[1024];
//						strcpy(filena_out_iter,argv[2]);
//						sprintf(iterna,"_%dth_iteration_output.tif",iter);
//						strcpy(&filena_out_iter[strlen(filena_out_iter)-4],iterna);
//
//						typedef itk::ImageFileWriter<ImageType> FileWriterType;
//						FileWriterType::Pointer writer2 = FileWriterType::New();
//						writer2->SetFileName(filena_out_iter);
//						writer2->SetInput(imout2);
//						try
//						{
//							writer2->Update();
//						}
//						catch (itk::ExceptionObject &e)
//						{
//							std::cout << e << std::endl;
//							return 1;
//						}
//
//				Vm2.clear();
//				//free(Cvotee);
//				printf("the average stick saliency is %f \n", means);
//				printf("%d-th iteration complete\n",iter);
//				printf("%ld seconds have elapsed since starting\n",time(NULL)-t1);
//				time_cal[iter+1]=time(NULL)-t1;
//
//			}
//
//			free(checked);
//			queuep.clear();
//
//			///////----------------------------------------end of tensor voting------------------------------------------------
//
//
//			//// Decompose the result for Densification.
//			//int rescounter=0;
//			//double sums=0;
//			//double means;
//			//double dividerv=0;
//			//double crit_0,crit_1,crit_2;
//			//for(int coz = window1;coz<pdepth-window1;coz++)
//			//	for(int coy = window; coy<pheight-window;coy+=1)
//			//		for(int cox = window; cox<pwidth-window;cox+=1)
//			//		{
//			//			if (C(coz,coy,cox)==1)
//			//			{
//			//				rescounter++;
//			//				continue;
//			//			}
//			//			votee_matrix=Vm2[rescounter].voting_matrix;
//			//			rescounter++;
//			//			rtvl_tensor<3> votee_tensor(votee_matrix);
//			//			crit_0=votee_tensor.saliency(0);
//			//			crit_1=votee_tensor.saliency(1);
//			//			crit_2=votee_tensor.saliency(2);
//			//			dividerv++;
//			//			sums=sums+crit_0;
//			//			printf("the stick saliency is %f \n", crit);
//			//			//if(crit_0>thre)
//			//			//if ((crit_0>crit_1)&&(crit_0>crit_2)&&(crit_0>thre))
//			//			if ((crit_0>crit_1)&&(crit_0>crit_2))
//			//			{
//			//				data Vres;
//			//				Vres.x=cox;
//			//				Vres.y=coy;
//			//				Vres.z=coz;
//			//				queuep.push_back(Vres);
//			//				C(coz,coy,cox)=1;
//			//			}
//
//			//		}
//			//		means=sums/dividerv;
//			//		printf("the average stick saliency is %f \n", means);
//
//
//
//			///////////////////////////////////////////////////////////////////////////
//			//
//			//end of tensor voting
//
//			//printf("About to print the data\n");
//			//ImageType::Pointer imout = ImageType::New();
//			//imout->SetRegions(im->GetLargestPossibleRegion());
//			//imout->Allocate();
//			//ImageType::IndexType output_index;
//			//for(int coz = window1;coz<pdepth-window1;coz++)
//			//	for(int coy = window; coy<pheight-window;coy+=1)// skip is used for sampling sometimes. We dont sample in Z direction
//			//		for(int cox = window; cox<pwidth-window;cox+=1)
//			//		{
//			//			if(C(coz,coy,cox)==1)
//			//			{
//			//				//int minimum_x=MAX(1,cox-window),minimum_y=MAX(1,coy-window),minimum_z=MAX(1,coz-window1);
//			//				//int maximum_x=MIN(rwidth/gridsize_x,cox+1-window),maximum_y=MIN(rlength/gridsize_y,coy+1-window),maximum_z=MIN(rdepth,coz+1-window1);
//
//			//				for(int cx = cox+1-window;cx<=cox+1-window;cx++)
//			//					for(int cy = coy+1-window;cy<=coy+1-window;cy++)
//			//					{
//			//						output_index[0] = xl+cx-1;
//			//						output_index[1] = yl+cy-1;
//			//						output_index[2] = zl+coz-window1;
//
//			//						fprintf(fp,"%d %d %d %d",zl+coz+1-window1,yl+cy,xl+cx,M(coz,coy,cox));
//			//						imout->SetPixel(output_index,255);
//			//						//#ifdef COMPUTE_L2
//			//						fprintf(fp," %lf",L(coz,coy,cox));
//			//						//#endif
//			//						fprintf(fp,"\n");
//			//					}
//
//			//			}
//			//		}
//			//typedef itk::ImageFileWriter<ImageType> FileWriterType;
//			//FileWriterType::Pointer writer = FileWriterType::New();
//			//writer->SetFileName(argv[2]);
//			//writer->SetInput(imout);
//			//try
//			//{
//			//	writer->Update();
//			//}
//			//catch (itk::ExceptionObject &e)
//			//{
//			//	std::cout << e << std::endl;
//			//	return 1;
//			//}
//
//			/*char distance_map_file[1024];
//			strcpy(distance_map_file,argv[2]);
//			strcpy(&distance_map_file[strlen(distance_map_file)-4],"_distance_map.tif");
//			char vector_map[1024];
//			strcpy(vector_map,argv[2]);
//			strcpy(&vector_map[strlen(vector_map)-4],"_vector_map.mhd");
//
//			typedef itk::Image<short int,3> DistanceImageType;
//			typedef itk::DanielssonDistanceMapImageFilter<ImageType,DistanceImageType> DistanceMapFilterType;
//			typedef DistanceMapFilterType::VectorImageType OffsetImageType;
//
//			DistanceMapFilterType::Pointer distfilter = DistanceMapFilterType::New();
//			distfilter->SetInput(imout);
//			distfilter->InputIsBinaryOn();
//			try
//			{
//				distfilter->Update();
//			}
//			catch (itk::ExceptionObject &e)
//			{
//				std::cout << e << std::endl;
//				return 1;
//			}
//
//			typedef itk::ImageFileWriter<DistanceImageType> DistFileWriter;
//			DistFileWriter::Pointer dwrite = DistFileWriter::New();
//			dwrite->SetFileName(distance_map_file);
//			dwrite->SetInput(distfilter->GetOutput());
//			try
//			{
//				dwrite->Update();
//			}
//			catch (itk::ExceptionObject &e)
//			{
//				std::cout << e << std::endl;
//				return 1;
//			}
//			
//			typedef itk::ImageFileWriter<OffsetImageType> OffsetFileWriter;
//			OffsetFileWriter::Pointer offwriter = OffsetFileWriter::New();
//			offwriter->SetFileName(vector_map);
//			offwriter->SetInput(distfilter->GetVectorDistanceMap());
//			try
//			{
//				offwriter->Update();
//			}
//			catch (itk::ExceptionObject &e)
//			{
//				std::cout << e << std::endl;
//				return 1;
//			}*/
//
//			//---comment these if reading segmented file directly
//			//free(matrix);
//			//free(lmatrix);
//			//free(checked);
//			//free(l1matrix);
//			//free(l2matrix);
//			//free(raster);
//			//fclose(fp);
//			for (int caltime=1;caltime<itert+1;caltime++)
//			{
//				printf("%d-th output took %ld seconds\n",caltime,time_cal[caltime]);
//			}
//			printf("%ld seconds \n",time(NULL)-t1);
//			return 0;
//}
