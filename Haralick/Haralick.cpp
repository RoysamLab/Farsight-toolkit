#include <iostream>
#include <math.h>
#include <time.h>
#include "Haralick.h"
#define SIGN(x,y) ((y)<0 ? -fabs(x) : fabs(x))
#ifdef _OPENMP
#include "omp.h"
#endif

using namespace std;


Haralick::Haralick(std::string imageFileName)
{
	this->imageFileName = imageFileName;		
    Ng=0;
	SetInputImage();
}


Haralick::Haralick(ImageType::Pointer inputImage)
{
	this->inputImage   = inputImage;
	IteratorType It( inputImage, inputImage->GetRequestedRegion());
	InputSizeType size = inputImage->GetRequestedRegion().GetSize();
	
	M=size[1];
	N=size[0];
	ImageType::IndexType requestedIndex =inputImage->GetRequestedRegion().GetIndex();

    f = (unsigned char **) malloc (M * sizeof(unsigned char *));
	for (int k = 0; k < M; k++)
		f[k] = (unsigned char *) malloc (N * sizeof(unsigned char));
	std::cout << "M = " << M << " N = " << N << std::endl;

    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
		ImageType::IndexType idx = It.GetIndex();            
        f[idx[1]][idx[0]]= It.Get(); 
	}
}


Haralick::Haralick()
{		
    Ng=0;
}
Haralick::~Haralick()
{

}



int Haralick:: ReadInputImage(std::string imageFileName)
{
	inputImage = ImageType::New();
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( imageFileName );
	try
	{
		reader->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught !" << std::endl; 
		std::cout << err << std::endl; 
		return -1;
	}
	inputImage = reader->GetOutput();

}

void Haralick::SetInputImage()
{

	ReadInputImage(imageFileName);
	IteratorType It( inputImage, inputImage->GetRequestedRegion());
	InputSizeType size = inputImage->GetRequestedRegion().GetSize();
	M=size[1];
	N=size[0];
	ImageType::IndexType requestedIndex =inputImage->GetRequestedRegion().GetIndex();

	f = (unsigned char **) malloc (M * sizeof(unsigned char *));
	for (int k = 0; k < M; k++)
		f[k] = (unsigned char *) malloc (N * sizeof(unsigned char));

	for (It.GoToBegin(); !It.IsAtEnd(); ++It)
	{

		ImageType::IndexType idx = It.GetIndex();
		f[idx[1]][idx[0]]= It.Get(); 
	}	

}


// this function calculates the maximum intensity value of image

int Haralick::max_value()
{
	int max=(f[0][0]);
	for(int i=0; i<M; i++)
	{
		for (int j=0; j<N; j++)
		{
			if((f[i][j])>max)
				max=f[i][j];
		}
	}

	return max;
}


// this method calculates minimum value of inensity value in image
int Haralick::min_value()
{
	int min=(f[0][0]);
	for(int i=0; i<M; i++)
	{
		for (int j=0; j<N; j++)
		{
			if((f[i][j])< min)
				min=f[i][j];
		}
	}
	if(min==0)
		return 1;
	else	
	    return min;
}
// the method sets the size of the cooccurence matrix 
void Haralick::set_Ng()
{
	Ng=(max_value()-min_value()+1);
}


//this method returns the no. of times a pattern repeates in the image (i.e. the (value1 value2) pattern in the image

//int Haralick::countrepetition(int value1,int value2,int dir[])
//{
//	int count=0;
//	int start_x=0,end_x=N-1,start_y=0,end_y=M-1;
//	///this was the original code for assymetric cooccurance matrix
//	if(dir[0]>=0 && dir[1]>=0)
//	{
//		start_x=0;
//		end_x=N-1-dir[0];
//		start_y=dir[1];
//		end_y=M-1;
//	}
//	else if(dir[0]>=0 && dir[1]<0)
//	{
//		start_x=0;
//		end_x=N-1-dir[0];
//		start_y=0;
//		end_y=M-1+dir[1];
//	}
//	else if(dir[0]<0 && dir[1]>=0)
//	{
//		start_x=-dir[0];
//		end_x=N-1;
//		start_y=dir[1];
//		end_y=M-1;
//	}
//	else 
//	{
//		start_x=-dir[0];
//		end_x=N-1;
//		start_y=0;
//		end_y=M-1+dir[1];
//	}
//	for(int i=start_y; i<=int(end_y); i++)
//	{
//		for (int j=start_x; j<=int(end_x); j++)
//		{
//			if((f[i][j]==value1) &&( f[(int)(i-dir[1])][(int)(j+dir[0])]==value2))
//				count++;
//		}
//	} 

	///////////////////////////////////////////////////////////
	/*if(dir[1]==0)
	{
	start_x=abs(dir[0]);
	end_x=N-1-abs(dir[0]);		

	}
	else if(dir[0]==0 )
	{	
	start_y=abs(dir[1]);
	end_y=M-1-abs(dir[1]);
	}
	else 
	{
	start_x=abs(dir[0]);
	end_x=N-1-abs(dir[0]);
	start_y=abs(dir[1]);
	end_y=M-1-abs(dir[1]);
	}

	for(int i=start_y; i<=int(end_y); i++)
	{
	for (int j=start_x; j<=int(end_x); j++)
	{
	if((f[i][j]==value1) &&( f[(int)(i-(dir[1]))][(int)(j+(dir[0]))]==value2))
	count=count+1;
	if((f[i][j]==value1) &&( f[(int)(i+(dir[1]))][(int)(j-(dir[0]))]==value2))
	count=count+1;

	}
	}     


	return count;
} */

/// now we need two intensity value
///for intensity1

/* this method gives the cooccurraance matrix given the dir[] vector
   user must be carefurl in deciding the dir[] vector
   the notation used is based on the cosine ans sine values of the directon angles
   cooccurance matrices in 8 directions 45* apart are calculated 
   for the current calculations of Haralick features the angle theta and theta+180  are considered same 
   i.e. to get a symmetric cooccurrence matrix the matrices in  theta and 180+theta direction are combined 
   for 0* the vectors are [1 0],[-1 0]
   for 90* the vectors are [0 1],[0 -1]
   for 45* the vectors are [1 1 ],[-1,-1]
   for 135* the vectors are [-1,1],[1,-1]
   please note  that the pixel separation in any dirextion(X or Y) is set 1 
   if the application requires more distance it ould be set accordingly
   */

int**  Haralick::calculateco_occurrence( int dir[])
{       
	    int** P2;
	    P2= ( int **) malloc (Ng * sizeof( int *));
			 for ( int k = 0; k < Ng; k++)
			 P2[k] = ( int *) malloc (Ng * sizeof( int));
		for(int i=0; i<Ng; i++)
		{
			for (int j=0; j<Ng; j++)
			{
				P2[i][j]=0; 
			}
		}
	    int start_x=0,end_x=N-1,start_y=0,end_y=M-1;
		if(dir[0]>=0 && dir[1]>=0)
		{
			start_x=0;
			end_x=N-1-dir[0];
			start_y=dir[1];
			end_y=M-1;
		}
		else if(dir[0]>=0 && dir[1]<0)
		{
			start_x=0;
			end_x=N-1-dir[0];
			start_y=0;
			end_y=M-1+dir[1];
		}
		else if(dir[0]<0 && dir[1]>=0)
		{
			start_x=-dir[0];
			end_x=N-1;
			start_y=dir[1];
			end_y=M-1;
		}
		else 
		{
			start_x=-dir[0];
			end_x=N-1;
			start_y=0;
			end_y=M-1+dir[1];
		}
		int min=min_value();
		int max=max_value();
		
		       
		for(int i=start_y; i<=int(end_y); i++)
		{
			for (int j=start_x; j<=int(end_x); j++)
			{
				if(f[i][j]>0)
				{
					int intensity1=f[i][j];
				    int intensity2=f[(int)(i-dir[1])][(int)(j+dir[0])];
				    P2[(intensity1-min)][(intensity2-min)]++;
				}
			}
		}        
		        
		return P2;

}

//this method calculates the cooccurance matrix for a given angle .As the cooccurrence matrix is symmetric 
//(used by haralick feature calculation) so the matrices coorsponding to angles theta and theta+180 are combined

void   Haralick::calculatecooccurrence( int** P ,int angle)
{

	int** P1;
	P1= ( int **) malloc (Ng * sizeof( int *));
	for ( int k = 0; k < Ng; k++)
		P1[k] = ( int *) malloc (Ng * sizeof( int));	
	int* dir1;int* dir2;
    

	dir1= ( int *) malloc (2 * sizeof( int ));
	dir2 = ( int *) malloc (2 * sizeof( int));	
    int min =min_value();

	if(angle==0)
	{
		dir1[0] = 1 ,dir1[1] = 0;                          //here cos0*=1 and sin0*=0
		dir2[0] = -1,dir2[1] = 0;
		
		P1=calculateco_occurrence( dir1);
		
		for(int i=0; i<Ng; i++)
		{
			for ( int j=0; j<Ng; j++)
			{
				P[i][j]=P1[i][j];
			}
		}
		P1=calculateco_occurrence( dir2);
		
		for(int i=0; i<Ng; i++)
		{
			for ( int j=0; j<Ng; j++)
			{
				P[i][j]=P[i][j]+P1[i][j];
			}
		}
       
	}
	
    else if(angle==45)
	{
		dir1[0]= 1,dir1[1]=1;                          // sines and cosines of directions are multiplied by sqrt(2) to acces the pixel at 45*
		dir2[0]= -1,dir2[1]= -1;
		
		P1=calculateco_occurrence( dir1);
		for(int i=0; i<Ng; i++)
		{
			for ( int j=0; j<Ng; j++)
			{
				P[i][j]=P1[i][j];
			}
		}
		P1=calculateco_occurrence( dir2);
		for(int i=0; i<Ng; i++)
		{
			for ( int j=0; j<Ng; j++)
			{
				P[i][j]=P[i][j]+P1[i][j];
			}
		}
        		
	}

	else if(angle==90)
	{
		dir1[0]=0 ,dir1[1]=1;                        
		dir2[0]= 0,dir2[1]= -1;
		P1=calculateco_occurrence( dir1);
		for(int i=0; i<Ng; i++)
		{
			for ( int j=0; j<Ng; j++)
			{
				P[i][j]=P1[i][j];
			}
		}
		P1=calculateco_occurrence( dir2);
		for(int i=0; i<Ng; i++)
		{
			for ( int j=0; j<Ng; j++)
			{
				P[i][j]=P[i][j]+P1[i][j];
			}
		}
		
	}
	else
	{
		dir1[0]= -1,dir1[1]=1;                                     // sines and cosines of directions are multiplied by sqrt(2)
		dir2[0]= 1,dir2[1]= -1;   
		P1=calculateco_occurrence( dir1);
		for(int i=0; i<Ng; i++)
		{
			for ( int j=0; j<Ng; j++)
			{
				P[i][j]=P1[i][j];
			}
		}
		P1=calculateco_occurrence( dir2);
		for(int i=0; i<Ng; i++)
		{
			for ( int j=0; j<Ng; j++)
			{
				P[i][j]=P[i][j]+P1[i][j];
			}
		}
	
	}
	/*for(int i=0; i<Ng; i++)
	{
		P[0][i] = 0;
		P[i][0] = 0;
	}*/
	}
/// this method calculate the haralick features 

void  Haralick::calculatefeatures(int** P )
{
	f1_ASM=0.0;f2_Contrast=0.0;f3_Correlation=0.0;f4_Variance=0.0;f5_IDM=0.0;f6_Sum_Avg=0.0;f7_Sum_Var=0.0;f8_Sum_Entropy=0.0;f9_Entropy=0.0;f10_diff_var=0.0;f11_difference_entropy=0.0;f12=0.0;f13=0.0;
	double R=0.0,mu=0.0,mu_x=0.0,mu_y=0.0,sigma_x=0.0,sigma_y=0.0;
	double HXY1=0.0,HXY2=0.0,HX=0.0,HY=0.0;
	
	int min=min_value();
	int max=max_value();
	int length1=2*(max-min)+1; ////length of p_xy1 array
	int length2=max-min+1;     ////// length of p_xy2 array

	double* p_x=(double *) malloc(Ng * sizeof(double));
	double* p_y=(double *) malloc(Ng * sizeof(double));
	double** p;
	p= ( double **) malloc (Ng * sizeof( double *));
	for ( int k = 0; k < Ng; k++)
		p[k] = ( double *) malloc (Ng * sizeof( double));

	for(int i=0;i<Ng;++i)
	{
		for(int j=0;j<Ng;++j)
		{
			p[i][j]=0;           

		}
	}   

	for(int i=0;i<Ng;++i)
	{
		p_x[i]=0.0;
		p_y[i]=0.0;
	}

	double* p_xy1 = (double *) malloc(length1* sizeof(double));
	double* p_xy2 = (double *) malloc(length2* sizeof(double));
	for(int i=0;i<length1;++i)
	{
		p_xy1[i] = 0.0;

	}
	for(int i=0;i<length2;++i)
	{
		p_xy2[i] = 0.0;

	}


	for(int i=0;i<Ng;++i)
	{
		for(int j=0;j<Ng;++j)
		{
			R+= P[i][j]; 
			

		}
	}
	
	for(int i=0;(i<Ng);++i)
	{
		for(int j=0;j<Ng;++j)
		{  
			p[i][j]= P[i][j]/R;      
			mu = mu + (p[i][j]*(i+min));
			p_y[i] += p[i][j];
		}
	}

	for(int i=0;(i<Ng);++i)
	{
		for(int j=0;j<Ng;++j)
		{  
			p_x[i]+= p[j][i];	
		}
	}

	

	for(int i=0;i<Ng;i++)
	{
		mu_x +=((i+min)*p_x[i]);
		mu_y +=((i+min)*p_y[i]);

		if(p_x[i]==0.0)
			HX -=0.0;
		else
			HX   -= (p_x[i]*log(p_x[i]));

		if(p_y[i]==0.0)
			HY-=0.0;
		else
			HY   -= (p_y[i]*log(p_y[i]));
	}

	for(int i=0;i<Ng;i++)
	{
		sigma_x +=(i-mu_x+min)*(i+min-mu_x)*(p_x[i]);
		sigma_y +=(i-mu_y+min)*(i+min-mu_y)*(p_y[i]);
	}



	//f1,f3_Correlation,f4_Variance,f5_IDM,f9_Entropy
	//double tmp_f1 = 0.0,  tmp_f5_IDM = 0.0, tmp_f9_Entropy = 0.0; 
	 double tmp_f3_Correlation = 0.0;
	
	//#pragma omp parallel for reduction(+:tmp_f1)
	for(int i=0;i<Ng;++i)
	{
		for(int j=0;j<Ng;++j)
		{
			f1_ASM+= (p[i][j]*p[i][j]);
			tmp_f3_Correlation += ((i+min)*(j+min)*p[i][j]);
			f5_IDM += (p[i][j]/(1+(i+min-j)*(i+min-j)));
			if(p[i][j]==0)
				f9_Entropy -=0.0;
			else
			f9_Entropy -=(p[i][j]*log(p[i][j]));
			f4_Variance += ((i-mu)*(i-mu)*p[i][j]);
		}
	}
	   f3_Correlation = (tmp_f3_Correlation-mu_x*mu_y)/(sigma_x*sigma_y);

	

	for(int n=0;n<=(max-min);++n)
	{  
		double l=0.0;
		for(int i=min;i<=max;++i)
		{
			for(int j=min; j<=max;++j)
			{
				if(abs(i-j)==n)
					l += p[i-min][j-min];
			}
		}
		f2_Contrast += (n*n*l);
	}           

	for(int i=0;i<Ng;++i)
	{
		for(int j=0;j<Ng;++j)
		{  
			if(p_x[i]==0.0 || p_y[j]==0)
			{   
				HXY1-=0.0;
				HXY2-=0.0;
			}
			else
			{
				HXY1 -=(p[i][j]*log(p_x[i]*p_y[j]));
				HXY2 -= (p_x[i]*p_y[j]*log(p_x[i]*p_y[j]));
			}
		}
	}


	for(int k=0;k<length1;++k)
	{  
		for(int i=min;i<=max;++i)
		{
			for(int j=min;(j<=max);++j)
			{
				if(i+j==(k+2*min)) 
					p_xy1[k]+=p[i-min][j-min];

			}
		}

	}

	for(int k=0;k<length2;++k)
	{  
		for(int i=min;i<=max;++i)
		{
			for(int j=min;j<=max;++j)
			{
				if(abs(i-j)==k)
					p_xy2[k]+=p[i-min][j-min];
			}

		}

	}          


	///f6_Sum_Avg,f7_Sum_Var,f8_Sum_Entropy

	for(int k=0;k<length1;++k)
	{
		f6_Sum_Avg += ((k+2*min)*p_xy1[k]);
		if(p_xy1[k]==0.0)
			f8_Sum_Entropy-=0.0;
		else
			f8_Sum_Entropy -=  ((p_xy1[k])*log(p_xy1[k]));
	}
	for(int k=0;k<length1;++k)
	{
		f7_Sum_Var += ((k-f6_Sum_Avg+2*min)*(k-f6_Sum_Avg+2*min)*(p_xy1[k]));
	}
	// f10_diff_var,f11_difference_entropy,ff
	double ff=0.0;
	for(int k=0;k<length2;++k)
	{
		ff += (k*p_xy2[k]);
		if(p_xy2[k]==0.0)
			f11_difference_entropy-=0.0;
		else
			f11_difference_entropy -= ((p_xy2[k])*log(p_xy2[k]));
	}
	for(int k=0;k<length2;++k)
	{

		f10_diff_var += (k-ff)*(k-ff)*p_xy2[k];
	}

	//f12,f13
	double m=HX;
	if(HX<HY)
		m = HY;

	f12=(f9_Entropy-HXY1)/m;
  	f13= sqrt(1-exp(-2.0*(HXY2-f9_Entropy)));
	////
	/////////EIGEN VECTORS------maximal largest eigen value


}


//////////////////////////////////////////////////////////
//lets mess up abit
double Haralick::f14_maxcorr (int** p,int  N)
 

/* Returns the Maximal Correlation Coefficient */
{
	//int i, j, k;
	double *px, *py, **Q;
	double *x, *iy;
	double tmp =0.0;
	double f1;
	px = pgm_vector (0, N);
	py = pgm_vector (0, N);
	Q = pgm_matrix (1, N + 1, 1, N + 1);
	x = pgm_vector (1, N);
	iy = pgm_vector (1, N);
	double** P;
	P= ( double **) malloc (N * sizeof( double *));
	for ( int k = 0; k < N; k++)
		P[k] = ( double *) malloc (N * sizeof( double));

	for(int i=0;i<N;++i)
	{
		for(int j=0;j<N;++j)
		{
			P[i][j]=0;           

		}
	}   

	for(int i=0;i<N;++i)
	{
		px[i]=0.0;
		py[i]=0.0;
	}

    double R=0.0;
	for(int i=0;i<N;++i)
	{
		for(int j=0;j<N;++j)
		{
			R+= p[i][j]; 		

		}
	}
	
	for(int i=0;(i<N);++i)
	{
		for(int j=0;j<N;++j)
		{  
			P[i][j]= p[i][j]/R;
			//cout<<P[i][j]<<" " ;
			
		}
	}



  /*
   * px[i] is the (i-1)th entry in the marginal probability matrix obtained
   * by summing the rows of p[i][j]
   */
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			px[i] += P[i][j];
			py[j] += P[i][j];
			//cout<<px[i]<<" " ;
		}
	}

  /* Find the Q matrix */
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			Q[i + 1][j + 1] = 0;
			for (int k = 0; k < N; ++k)
			{
				if(px[i]!=0 && py[k] != 0)
					Q[i + 1][j + 1] += P[i][k] * P[j][k] / (px[i] * py[k]);
				else
					Q[i + 1][j + 1] += 0;
			}

			//cout<<Q[i+1][j+1]<<" " ;
		}
	}

  /* Balance the matrix */
	mkbalanced (Q, N);
  /* Reduction to Hessenberg Form */
	reduction (Q, N);
  /* Finding eigenvalue for nonsymetric matrix using QR algorithm */
	if (!hessenberg (Q,  x, iy,N))
	{
		for (int i=1; i<=N+1; i++)
			free(Q[i]+1);
		free(Q+1);
		free((char *)px);
		free((char *)py);
		free((x+1));
		free((iy+1));
		return 0.0;
	  
	}
   simplesrt(N,x); 
  /* Returns the sqrt of the second largest eigenvalue of Q */
	for (int i = 2, tmp = x[1]; i <= N; ++i)
	{
		tmp = (tmp > x[i]) ? tmp : x[i];
		//cout<<tmp<< " " ;
	}
	f1 = sqrt(x[N - 1]);
	for (int i=1; i<=N+1; i++) 
		free(Q[i]+1);
	free(Q+1);
	free((char *)px);
	free((char *)py); 
	free((x+1)); 
	free((iy+1));
	  return f1;
}

double *Haralick::pgm_vector (int nl, int nh)
 
{
  double *v;
  int    i;

  v = (double *) malloc ((unsigned) (nh - nl + 1) * sizeof (double));
  if (!v)
    fprintf (stderr, "memory allocation failure (pgm_vector) "), exit (1);

  for (i=0; i<=(nh-nl); i++) v[i]=0;
  return v - nl;
}


double **Haralick::pgm_matrix (int nrl,int  nrh, int ncl, int nch)


/* Allocates a double matrix with range [nrl..nrh][ncl..nch] */
{
  int i;
  double **m;

  /* allocate pointers to rows */
  m = (double **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (double *));
  if (!m)
    fprintf (stderr, "memory allocation failure (pgm_matrix 1) "), exit (1);
  m -= ncl;

  /* allocate rows and set pointers to them */
  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (double *) malloc ((unsigned) (nch - ncl + 1) * sizeof (double));
    if (!m[i])
      fprintf (stderr, "memory allocation failure (pgm_matrix 2) "), 
	exit (2);
    m[i] -= ncl;
  }
  /* return pointer to array of pointers to rows */
  return m;
}

void Haralick::results (double *Tp,char *c, double *a)

{
  int i;
  double max, min;
  max = a[0];
  min = a[0];
/*  DOT;
  fprintf (stdout, "%s", c);
*/  for (i = 0; i < 4; ++i, *Tp++)
    {	
    if (a[i] <= min)
	min = a[i];
    if (a[i] > max)
	max = a[i];
  /*  fprintf (stdout, "% 1.3e ", a[i]); */
    *Tp = a[i];
    }	
/*  fprintf (stdout, "% 1.3e  % 1.3e\n", 
    (a[0] + a[1] + a[2] + a[3]) / 4,max-min); */
  *Tp = (a[0] + a[1] + a[2] + a[3]) / 4;
  *Tp++;
  *Tp = max - min;
 
  	
}

void Haralick::simplesrt ( int n,  double arr[])

{
  int i, j;
  double a;

  for (j = 2; j <= n; j++)
  {
    a = arr[j];
    i = j - 1;
    while (i > 0 && arr[i] > a)
    {
      arr[i + 1] = arr[i];
      i--;
    }
    arr[i + 1] = a;
  }
}

void Haralick::mkbalanced (  double **a,  int n)

{
  int last, j, i;
  double s, r, g, f, c, sqrdx;
  double RADIX = 2.0;

  sqrdx = RADIX * RADIX;
  last = 0;
  while (last == 0)
  {
    last = 1;
    for (i = 1; i <= n; i++)
    {
      r = c = 0.0;
      for (j = 1; j <= n; j++)
	if (j != i)
	{
	  c += fabs (a[j][i]);
	  r += fabs (a[i][j]);
	}
      if (c && r)
      {
	g = r / RADIX;
	f = 1.0;
	s = c + r;
	while (c < g)
	{
	  f *= RADIX;
	  c *= sqrdx;
	}
	g = r * RADIX;
	while (c > g)
	{
	  f /= RADIX;
	  c /= sqrdx;
	}
	if ((c + r) / f < 0.95 * s)
	{
	  last = 0;
	  g = 1.0 / f;
	  for (j = 1; j <= n; j++)
	    a[i][j] *= g;
	  for (j = 1; j <= n; j++)
	    a[j][i] *= f;
	}
      }
    }
  }
}


void Haralick::reduction ( double **a,  int n)

{
  int m, j, i;
  double y, x;

  for (m = 2; m < n; m++)
  {
    x = 0.0;
    i = m;
    for (j = m; j <= n; j++)
    {
      if (fabs (a[j][m - 1]) > fabs (x))
      {
		  x = a[j][m - 1];
		  i = j;
      }
    }
    if (i != m)
    {
      for (j = m - 1; j <= n; j++)
		  SWAP (a[i][j], a[m][j]);
	  for (j = 1; j <= n; j++)
		  SWAP (a[j][i], a[j][m]) ;
	  a[j][i] = a[j][i];
    }
    if (x)
    {
      for (i = m + 1; i <= n; i++)
      {
	if (y = a[i][m - 1])
	{
	  y /= x;
	  a[i][m - 1] = y;
	  for (j = m; j <= n; j++)
	    a[i][j] -= y * a[m][j];
	  for (j = 1; j <= n; j++)
	    a[j][m] += y * a[j][i];
	}
      }
    }
  }
}

int Haralick:: hessenberg ( double **a,double wr[],double wi[],  int n)
 

{
  int nn, m, l, k, j, its, i, mmin;
  double z, y, x, w, v, u, t, s, r, q, p, anorm;

  anorm = fabs (a[1][1]);
  for (i = 2; i <= n; i++)
    for (j = (i - 1); j <= n; j++)
      anorm += fabs (a[i][j]);
  nn = n;
  t = 0.0;
  while (nn >= 1)
  {
    its = 0;
    do
    {
      for (l = nn; l >= 2; l--)
      {
	s = fabs (a[l - 1][l - 1]) + fabs (a[l][l]);
	if (s == 0.0)
	  s = anorm;
	if ((double) (fabs (a[l][l - 1]) + s) == s)
	  break;
      }
      x = a[nn][nn];
      if (l == nn)
      {
	wr[nn] = x + t;
	wi[nn--] = 0.0;
      }
      else
      {
	y = a[nn - 1][nn - 1];
	w = a[nn][nn - 1] * a[nn - 1][nn];
	if (l == (nn - 1))
	{
	  p = 0.5 * (y - x);
	  q = p * p + w;
	  z = sqrt (fabs (q));
	  x += t;
	  if (q >= 0.0)
	  {
	    z = p + SIGN (z, p); 
	    wr[nn - 1] = wr[nn] = x + z;
	    if (z)
	      wr[nn] = x - w / z;
	    wi[nn - 1] = wi[nn] = 0.0;
	  }
	  else
	  {
	    wr[nn - 1] = wr[nn] = x + p;
	    wi[nn - 1] = -(wi[nn] = z);
	  }
	  nn -= 2;
	}
	else
	{
	  if (its == 30)
	    {
/*	    fprintf (stderr, 
  "Too many iterations to required to find %s\ngiving up\n", F14);  */
	     return 0; /*exit (1);*/
	     }			
	  if (its == 10 || its == 20)
	  {
	    t += x;
	    for (i = 1; i <= nn; i++)
	      a[i][i] -= x;
	    s = fabs (a[nn][nn - 1]) + fabs (a[nn - 1][nn - 2]);
	    y = x = 0.75 * s;
	    w = -0.4375 * s * s;
	  }
	  ++its;
	  for (m = (nn - 2); m >= l; m--)
	  {
	    z = a[m][m];
	    r = x - z;
	    s = y - z;
	    p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
	    q = a[m + 1][m + 1] - z - r - s;
	    r = a[m + 2][m + 1];
	    s = fabs (p) + fabs (q) + fabs (r);
	    p /= s;
	    q /= s;
	    r /= s;
	    if (m == l)
	      break;
	    u = fabs (a[m][m - 1]) * (fabs (q) + fabs (r));
	    v = fabs (p) * (fabs (a[m - 1][m - 1]) + 
			    fabs (z) + fabs (a[m + 1][m + 1]));
	    if ((double) (u + v) == v)
	      break;
	  }
	  for (i = m + 2; i <= nn; i++)
	  {
	    a[i][i - 2] = 0.0;
	    if (i != (m + 2))
	      a[i][i - 3] = 0.0;
	  }
	  for (k = m; k <= nn - 1; k++)
	  {
	    if (k != m)
	    {
	      p = a[k][k - 1];
	      q = a[k + 1][k - 1];
	      r = 0.0;
	      if (k != (nn - 1))
		r = a[k + 2][k - 1];
	      if (x = fabs (p) + fabs (q) + fabs (r))
	      {
		p /= x;
		q /= x;
		r /= x;
	      }
	    }
	    if (s = SIGN (sqrt (p * p + q * q + r * r), p)) 
	    {
	      if (k == m)
	      {
		if (l != m)
		  a[k][k - 1] = -a[k][k - 1];
	      }
	      else
		a[k][k - 1] = -s * x;
	      p += s;
	      x = p / s;
	      y = q / s;
	      z = r / s;
	      q /= p;
	      r /= p;
	      for (j = k; j <= nn; j++)
	      {
		p = a[k][j] + q * a[k + 1][j];
		if (k != (nn - 1))
		{
		  p += r * a[k + 2][j];
		  a[k + 2][j] -= p * z;
		}
		a[k + 1][j] -= p * y;
		a[k][j] -= p * x;
	      }
	      mmin = nn < k + 3 ? nn : k + 3;
	      for (i = l; i <= mmin; i++)
	      {
		p = x * a[i][k] + y * a[i][k + 1];
		if (k != (nn - 1))
		{
		  p += z * a[i][k + 2];
		  a[i][k + 2] -= p * r;
		}
		a[i][k + 1] -= p * q;
		a[i][k] -= p;
	      }
	    }
	  }
	}
      }
    } while (l < nn - 1);
  }
return 1;
}

void  Haralick::SWAP(double a,double b)
{
	double y;
	y=a;
	a=b;
	b=y;
}
double * Haralick::GetHaralik()
{
	
	set_Ng();
	double * features = (double *) malloc(14 * sizeof(double)); 
	for(int k=0;k<14;++k)
		features[k]=0.0;

	double* f1= (double *) malloc(14 * sizeof(double));
	int** P;
	P= ( int **) malloc ( Ng * sizeof( int *));
	for ( int k = 0; k < Ng; k++)
		P[k] = ( int *) malloc (Ng * sizeof( int));
	for(int i=0;i<Ng;i++)
	{
		for(int j=0;j<Ng;j++)
		{
			P[i][j]=0;
		}
	}
    
	//#pragma omp parallel for
	for(int angle=0;angle<=135;angle+=45)
	{
		calculatecooccurrence(P,angle);
		calculatefeatures(P);
		f1[0]=f1_ASM;f1[1]=f2_Contrast;f1[2]=f3_Correlation;f1[3]=f4_Variance;f1[4]=f5_IDM;f1[5]=f6_Sum_Avg;f1[6]=f7_Sum_Var;
		f1[7]=f8_Sum_Entropy;f1[8]=f9_Entropy;f1[9]=f10_diff_var;f1[10]=f11_difference_entropy;f1[11]=f12;f1[12]=f13;
		f1[13]=f14_maxcorr (P, Ng);

		for(int k=0;k<14;++k)
		{
			features[k]+=f1[k]/4.0;
		}
	}
	
	return features;
}