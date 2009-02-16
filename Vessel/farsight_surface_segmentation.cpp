//#include <stdio.h>
#include "farsight_surface_segmentation.h"

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

//variables
// int pdepth,pheight,pwidth;
// PixelType ** lmatrix, **matrix;
// int window, int window1;
// double alpha1;
#define DEBUG_ printf



#define CFH(a,b) (((unsigned int)a)<<(b*4))


void farsight_surface_segmentation::set_image(std::string const& image_path,std::string const& image_name)
{
	image_filename = image_path+image_name;
}

void farsight_surface_segmentation::set_parameters(std::string const& xml_filename)
{
	string copy = xml_filename;
	params_surface parsed = parseXML_surface((char*)(copy.c_str()),parsed); // file overwrites what is passed to the function
	skip = parsed.skip;
	window = parsed.window;
	window1 = parsed.window1;
	alpha1 = parsed.alpha1;
	epsilonw = parsed.epsilonw;
	debug = parsed.debug;
}

return_data farsight_surface_segmentation::function(int coz, int coy, int cox)
{
	//printf("Entered function %d %d %d pwidth %d pimsize %d\n",coz, coy, cox,pwidth,pimsize);
	double t[3];
	double best_t1= -1;
	double best_t0 = -1;
	//the main for loop
	double L1 = -5000000;
	double L;
	double amax,bmax,cmax;
	int pcount =0;
	int pccn =0;
	//printf("starting loop for null\n");
		for(int wz=-window1;wz<=window1;wz++)
			for(int wy=-window;wy<=window;wy++)
				for(int wx=-window;wx<=window;wx++)
				{
				//	printf("pccn %d coz+wz %d coy+wy %d cox+wx %d\n",pccn,coz+wz,coy+wy,cox+wx);
					null[pccn++]=P(coz+wz,coy+wy,cox+wx);
				}
//	printf("Ending loop for null\n");
	//defined in stdlib to sort sequential arrays
	qsort(null,pccn,sizeof(unsigned char),compare);
	// find the median of the null hypothesis
	/*double sum_sort =0;
	int sum_sort_num =0;
	for(int counter_s=alpha1*pccn; counter_s<pccn-alpha1*pccn;counter_s++,sum_sort_num++)
	sum_sort+=null[counter_s];
	t[2]=sum_sort/sum_sort_num;
	*/
	if(pccn%2==0)
		t[2]=(null[pccn/2]+null[pccn/2-1])/2.0;
	else
		t[2]=null[pccn/2];

	// a small threshold fixed to accelerate the program. Anything with a median of less than 2 is ignored

	if(t[2]<prune)
	{
		return_data ret;
		ret.lvalue = -10;
		ret.M = -10;
	//	printf("I came here before return ret\n");
		return ret;
		//continue;
	}
	double a,b,c;
	//printf("I did come here inside the function starting alternate hypothesis\n");
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
					if((wz*a+wy*b+wx*c)<= epsilonw && (wz*a+wy*b+wx*c)>=-epsilonw)
						fore[pccf++]=P(coz+wz,coy+wy,cox+wx);
					else
						back[pccb++]=P(coz+wz,coy+wy,cox+wx);
				}

				//this is an inefficient way to find the median. 
				//I sort the arrays and find the middle element
				qsort(fore,pccf,sizeof(unsigned char),compare);
				qsort(back,pccb,sizeof(unsigned char),compare);

				//find the median of back
				if(pccb%2==0)
					t[0]=(back[pccb/2]+back[pccb/2-1])/2.0;
				else
					t[0]=back[pccb/2];
				//find the median of fore
				if(pccf%2==0)
					t[1]=(fore[pccf/2]+fore[pccf/2-1])/2.0;
				else
					t[1]=fore[pccf/2];

				//compare the difference of two lambdas to a lambda_threshold =5
				double likelihood=0;
				if(t[1]-t[0]>=lambda_lower_threshold)
				{
					L=t[1]-t[0];// we dont actually compute the L because its not needed. This would render the normal vectors useless but we dont need them anyway.

															if(t[2]==0)
										{
											likelihood == 10000;
										}
										else
										{
											likelihood += t[2]-t[1] +t[1]*log(t[1]/t[2]);
											likelihood += t[2]-t[0];
											if(t[0]!=0)
												likelihood+=t[0]*log(t[0]/t[2]);

										}

					//double sum_sort =0;
					//int sum_sort_num =0;
					//for(int counter_s=int(alpha1*pccf); counter_s<int(pccf-alpha1*pccf);counter_s++)
					//	sum_sort+=fore[counter_s];

					//if(t[1]>0)
					//	likelihood+=-pccf*(1-2*alpha1)*t[1]+sum_sort*log(t[1]);
					//if(t[2]>0)
					//	likelihood-=(-pccf*(1-2*alpha1)*t[2]+sum_sort*log(t[2]));
					//sum_sort =0;
					//sum_sort_num =0;
					//for(int counter_s=alpha1*pccb; counter_s<pccb-alpha1*pccb;counter_s++)
					//	sum_sort+=back[counter_s];
					//if(t[0]>0)
					//	likelihood +=-pccb*(1-2*alpha1)*t[0]+sum_sort*log(t[0]);
					//if(t[2]>0)
					//	likelihood -=(-pccb*(1-2*alpha1)*t[2]+sum_sort*log(t[2]));



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

	return_data ret;
	ret.M = int(best_t1-best_t0+0.5);
	ret.lvalue = L1;
//	printf("Leaving function\n");
	return ret;

}

bool farsight_surface_segmentation::run()
{
	time_t t1,t2;
	t1 = time(NULL);
	unsigned char * raster;
	
	int rwidth,rlength,rdepth;
	// pic file name

	
	FILE * fpi = fopen(image_filename.c_str(),"rb");
	if(fpi == NULL)
	{
		std::cerr<<"I tried this file "<<image_filename.c_str()<<std::endl;
		std::cerr<<"Couldn't open the file. Please check the filename and try again"<<std::endl;
		return false;
	}

	char *temp=new char[77];
	fread(temp,sizeof(char),76,fpi);
	//printf("%d %d %d\n",CFH(temp[0],0)+CFH(temp[1],2),CFH(temp[2],0)+CFH(temp[3],2),CFH(temp[4],0)+CFH(temp[5],2));
	/*for(int counter =0; counter<10; counter++)
	{
	printf("%d %d\n",counter,temp[counter]);
	}*/
	rwidth = CFH(temp[0],0)+CFH(temp[1],2);
	rlength = CFH(temp[2],0)+CFH(temp[3],2);
	rdepth = CFH(temp[4],0)+CFH(temp[5],2);
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
	//printf("I'm going to read the memory block now rdepth = %d rheight = %d rwidth = %d\n", rdepth, rlength, rwidth);
	int number=0;
	int sum = npixels*rdepth;
	//do{
	number = fread(raster+npixels*rdepth-sum,sizeof(unsigned char),sum-number,fpi);
	sum = sum-number;
	printf("done reading\n");
	if(feof(fpi))
	{
		printf("got feof\n");
	}
	clearerr(fpi);
	printf("number of elements read = %d\n",number);
	//}while(number>0);
	fclose(fpi);
	printf("values %d %d %d\n",rdepth, rwidth, rlength);



	/*
	char filenameoutbuff[1024];
	sprintf(filenameoutbuff,"%s.npts",argv[1]);
	FILE * fp = fopen(filenameoutbuff,"w");*/
	
	//double *lmatrix;
	//int * matrix;
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
		-1,    -1,    -1,
			-1,     0,    -1,
			-1,     1,    -1,
			0,    -1,    -1,
			0,     0,    -1,
			0,     1,    -1,
			1,    -1,    -1,
			1,     0,    -1,
			1,     1,    -1,
			-1,    -1,     0,
			-1,     0,     0,
			-1,     1,     0,
			0,    -1,     0,
			0,     1,     0,
			1,    -1,     0,
			1,     0,     0,
			1,     1,     0,
			-1,    -1,     1,
			-1,     0,     1,
			-1,     1,     1,
			0,    -1,     1,
			0,     0,     1,
			0,     1,     1,
			1,    -1,     1,
			1,     0,     1,
			1,     1,     1
	};
	// the following control the limits in the input image that need to be segmented. index starts from 0
	//rlength = 512;
	//rwidth = 512;
	//Srdepth=26;
	//for(int gridx = 0; gridx <gridsize_x; gridx++)
	//{
	//	for(int gridy = 0; gridy <gridsize_y;gridy++)
	//	{
	//s(:,900:1024,250:350
	//int zl = 0,zu=rdepth-1,yl=899,yu=1023,xl=249,xu=349;
	int zl = 0,zu=rdepth-1,yl=gridy*(rlength/gridsize_y),yu=(gridy+1)*(rlength/gridsize_y)-1,xl=gridx*(rwidth/gridsize_x),xu=(gridx+1)*(rwidth/gridsize_x)-1;
	//	printf("starting for grid location %d %d\n",gridy+1,gridx+1);
	//set_limits(zl,zu,yl,yu,xl,xu,argc,argv);
	printf("zl zu yl yu xl xu = %d %d %d %d %d %d\n",zl,zu,yl,yu,xl,xu);

	//int zl = 0,zu=57,yl=0,yu=511,xl=0,xu=511;

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
	// just zero padded as of now
	for (int coz=window1;coz<pdepth-window1;coz++)
		for(int coy=window;coy<pheight-window;coy++)
			for(int cox=window;cox<pwidth-window;cox++)
			{
				/*if(MAT(zl+coz-window1,yl+coy-window,xl+cox-window)>0)
				{
				printf("mat was >0");
				}*/
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
			printf("Did I come till epsilon?\n");
			// continue adding more variables
			//						double alpha1 = 0.3;
			double epsilon = 1e-2;

			//Filename used for writing the output

			//unsigned int pc = 0;

			double scaling_matrix[] = {window1, window, window};
			// 21 directions are defined and they are reflected to give 42 directions
			double planes_21[21][3]= { 1, 0, -1,
				1, 0, 1,
				1, -1, 0,
				1, 1, 0,
				0, 1, 1,
				0, -1, 1,
				1, 0, 0,
				0, 1, 0,
				0, 0, 1,
				0.5, 0.5, 1,
				0.5, -0.5, 1,
				-0.5, -0.5, 1,
				-0.5, 0.5, 1,
				1, 0.5, 0.5,
				1, 0.5, -0.5,
				1, -0.5, 0.5,
				1, -0.5, -0.5,
				0.5, 1, 0.5,
				-0.5, 1, 0.5,
				0.5, 1, -0.5,
				-0.5, 1, -0.5};

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
			epsilonw=1.5;
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
						int pcount =0;
						int pccn =0;

						for(int wz=-window1;wz<=window1;wz++)
							for(int wy=-window;wy<=window;wy++)
								for(int wx=-window;wx<=window;wx++)
									null[pccn++]=P(coz+wz,coy+wy,cox+wx);

						//defined in stdlib to sort sequential arrays
						qsort(null,pccn,sizeof(unsigned char),compare);
						// find the median of the null hypothesis
						/*double sum_sort =0;
						int sum_sort_num =0;
						for(int counter_s=alpha1*pccn; counter_s<pccn-alpha1*pccn;counter_s++,sum_sort_num++)
						sum_sort+=null[counter_s];
						t[2]=sum_sort/sum_sort_num;
						*/
						if(pccn%2==0)
							t[2]=(null[pccn/2]+null[pccn/2-1])/2.0;
						else
							t[2]=null[pccn/2];

						// a small threshold fixed to accelerate the program. Anything with a median of less than 2 is ignored

						if(t[2]<prune)
							continue;
						double a,b,c;
						//	printf("I did come here\n");
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
										if((wz*a+wy*b+wx*c)<= epsilonw && (wz*a+wy*b+wx*c)>=-epsilonw)
											fore[pccf++]=P(coz+wz,coy+wy,cox+wx);
										else
											back[pccb++]=P(coz+wz,coy+wy,cox+wx);
									}

									//this is an inefficient way to find the median. 
									//I sort the arrays and find the middle element
									qsort(fore,pccf,sizeof(unsigned char),compare);
									qsort(back,pccb,sizeof(unsigned char),compare);

									//find the median of back
									if(pccb%2==0)
										t[0]=(back[pccb/2]+back[pccb/2-1])/2.0;
									else
										t[0]=back[pccb/2];
									//find the median of fore
									if(pccf%2==0)
										t[1]=(fore[pccf/2]+fore[pccf/2-1])/2.0;
									else
										t[1]=fore[pccf/2];

									//compare the difference of two lambdas to a lambda_threshold =5
									likelihood=0;
									if(t[1]-t[0]>=lambda_lower_threshold)
									{
										L=t[1]-t[0];// we dont actually compute the L because its not needed. This would render the normal vectors useless but we dont need them anyway.

										if(t[2]==0)
										{
											likelihood == 10000;
										}
										else
										{
											likelihood += t[2]-t[1] +t[1]*log(t[1]/t[2]);
											likelihood += t[2]-t[0];
											if(t[0]!=0)
												likelihood+=t[0]*log(t[0]/t[2]);

										}
										//double sum_sort =0;
										//int sum_sort_num =0;
										//for(int counter_s=alpha1*pccf; counter_s<pccf-alpha1*pccf;counter_s++,sum_sort_num++)
										//	sum_sort+=fore[counter_s];

										//if(t[1]>0)
										//	likelihood+=-pccf*(1-2*alpha1)*t[1]+sum_sort*log(t[1]);
										//if(t[2]>0)
										//	likelihood-=(-pccf*(1-2*alpha1)*t[2]+sum_sort*log(t[2]));
										//sum_sort =0;
										//sum_sort_num =0;
										//for(int counter_s=alpha1*pccb; counter_s<pccb-alpha1*pccb;counter_s++,sum_sort_num++)
										//	sum_sort+=back[counter_s];
										//if(t[0]>0)
										//	likelihood +=-pccb*(1-2*alpha1)*t[0]+sum_sort*log(t[0]);
										//if(t[2]>0)
										//	likelihood -=(-pccb*(1-2*alpha1)*t[2]+sum_sort*log(t[2]));


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
							double L2 = 0;

						}
						else
						{
							M(coz,coy,cox)=-1;

							//#ifdef COMPUTE_L2
							//										L(coz,coy,cox)=L1;
							//#endif
						}
					}
					printf("%d coz\n",coz);
					printf("%d seconds have elapsed since starting\n",time(NULL)-t1);
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
				if(pc == queuep.size())
					break;
				data seedp = queuep[pc];
				assert(M(seedp.z,seedp.y,seedp.x)>0);
				int l1 = M(seedp.z,seedp.y,seedp.x);

				double l2 = L(seedp.z,seedp.y,seedp.x);

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
					return_data return_value;
					if(M(z1,y1,x1)==0)
					{
						return_value = function(z1,y1,x1);
						M(z1,y1,x1) = return_value.M;
						L(z1,y1,x1) = return_value.lvalue;
						if(L(z1,y1,x1) <0)
						{
							M(z1,y1,x1) = -1;// null this if L <0 
						}
					}
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

					free(checked);
					free(l1matrix);
					free(l2matrix);
					//	}
					//	}

					free(raster);
				
					printf("%d seconds \n",time(NULL)-t1);
					return true;
}

void farsight_surface_segmentation::read_xml(std::string const& file_path,std::string const& filename)
{
	free(matrix);
	free(lmatrix);
	string output_filename = file_path+filename;
	params_surface par;
	bool success = parseXML_output((char*)output_filename.c_str(),matrix,lmatrix,par);
	imptr = ImageType::New();
	ImageType::SizeType size;
	size[0] = par.sizez;
	size[1] = par.sizey;
	size[2] = par.sizex;
	ImageType::IndexType index;
	index[0] = 0;
	index[1] = 0;
	index[2] = 0;
	ImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(index);
	imptr->SetRegions(region);
	imptr->Allocate();
	
	for(int coz = window1;coz<pdepth-window1;coz++){
		index[0]= coz-window1;
		for(int coy = window; coy<pheight-window;coy+=1){
			index[1] = coy-window;
			for(int cox = window; cox<pwidth-window;cox+=1){
				index[2] = cox - window;
				if(M(coz,coy,cox)>0){
					imptr->SetPixel(index,255);
				}
			}
		}
	}

}
farsight_surface_segmentation::ImageConstPtrType farsight_surface_segmentation::display()
{
	return imptr.GetPointer();
}
void farsight_surface_segmentation::write_xml(std::string const& file_path,std::string const& xml_filename)
 {

	 string output_filename = file_path+xml_filename;
	 FILE * fp = fopen(output_filename.c_str(),"w");

	printf("About to print the data\n");
	fprintf(fp,"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
	fprintf(fp,"<vessel program= \"surface_segmentation\" version = \"1.0\" >\n");
	fprintf(fp,"<parameters sizex=\"%d\" sizey=\"%d\" sizez=\"%d\" window_z=\"%\" window_xy=\"%d\" foreground_width=\"%d\" do_debug=\"%d\" alpha=\"%d\"\\>\n",pwidth-2*window,pheight-2*window,pdepth-2*window1,window1,window,epsilonw,debug,alpha1);
	fprintf(fp,"<segmentation>\n");
	for(int coz = window1;coz<pdepth-window1;coz++)
		for(int coy = window; coy<pheight-window;coy+=1)// skip is used for sampling. We dont sample in Z direction
				for(int cox = window; cox<pwidth-window;cox+=1)
					if(M(coz,coy,cox)>0)
					{
						fprintf(fp,"<point z = \"%d\" y = \"%d\" x = \"%d\" d=\"%d\" l=\"%lf\"\\>\n",coz+1-window1,coy-window+1,cox-window+1,M(coz,coy,cox),L(coz,coy,cox));
					}
	fprintf(fp,"<\\segmentation>\n");
	fclose(fp);
}





