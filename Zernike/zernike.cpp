#include "zernike.h"

///////////////////////
/*THE CLASS ZERNIKE RETURNS THE ZERNIKE MOMENTS FOR GRAY SCALE IMAGES
Now,zernike moments are orthogonal on the unit circle(r==1) */
////////////////////////////
namespace ftk{
zernike::zernike()
{
	M =0;
    N =0;

}
zernike::zernike(std::string imageFileName,int orderofmoments)
{
	this->imageFileName   = imageFileName;
	this ->orderofmoments = orderofmoments;
	M =0;
    N =0;
	SetInputImage();

}
zernike::zernike(ImageType::Pointer inputImage,int orderofmoments)
{
	this->inputImage   = inputImage;
	this ->orderofmoments = orderofmoments;
	
	IteratorType It( inputImage, inputImage->GetRequestedRegion());
	InputSizeType size = inputImage->GetRequestedRegion().GetSize();
	
	M=size[1];
	N=size[0];
	ImageType::IndexType requestedIndex =inputImage->GetRequestedRegion().GetIndex();

    f = (unsigned char **) malloc (M * sizeof(unsigned char *));
	for (int k = 0; k < M; k++)
		f[k] = (unsigned char *) malloc (N * sizeof(unsigned char));
	//std::cout << "M = " << M << " N = " << N << std::endl;

    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
		ImageType::IndexType idx = It.GetIndex();            
        f[idx[1]][idx[0]]= It.Get(); 
	}
}




zernike::~zernike()
{

}

/* this method is used to set the input image parameters using  Image Region Iterators
as we need image dimensions (M and N)  for the calculations of Zernike moments and need to access the image pixel data */

int zernike:: ReadInputImage(std::string imageFileName)
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

void zernike::SetInputImage()
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
	//std::cout << "M = " << M << " N = " << N << std::endl;

    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
		ImageType::IndexType idx = It.GetIndex();            
        f[idx[1]][idx[0]]= It.Get(); 
		/*if(int(It.Get()>0))
			f[idx[1]][idx[0]]=255;*/
	}  
}

/* this methog gives the values of X coordinates in the  image*/
double* zernike::getX()
{
	double* X = (double *) malloc(N * sizeof(double));
    for(int i=0;i<=N-1;++i)
    {
		
		X[i] = i;
		
	}
	return X;
}
/* this methog gives the values of Y coordinates in the mapped image*/
double* zernike::getY()
{
	double* Y = (double *) malloc(M * sizeof(double));     
    for(int i=0;i<=M-1;++i)
    {
	   Y[i] = i;
    }
    return Y;
}

/*this method returns the (0,0) order geometrical moment
which is equal to sum(sum(f(x,y));*/
double zernike::G00()
{
	double temp =0.0;
	for(int i=0;i<M;++i)
	{
		for(int j=0;j<N;++j)
		{
			temp+=f[i][j];
		}
    }
	return (temp);
	
}
/*this method returns the (1,0) order geometrical moment which is needed to calculate centre of floroscence for the image
( i.e. X component of  centre of mass in image is is G10()/G00()); */

double zernike::COF_X()
{
	double Moments=0.0;
	double* X =getX();              
	for(int i=0;i<M;++i)
	{
		for(int j=0;j<N;++j)
		{
			Moments+=((f[i][j])*X[j]);
		}
	}
	free(X);
	Moments=Moments/G00();
	return Moments;
}
/*this method returns the (0,1) order geometrical moment which is needed to calculate the  centre of fluoroscence for the image
( i.e. Y component of  centre of mass in image is is G01()/G00()); */


double zernike::COF_Y()
{
	double Moments=0.0;
	double* Y = getY();          
	for(int i=0;i<M;++i)
	{
		for(int j=0;j<N;++j)
		{
			Moments+=((f[i][j])*Y[i]);
		}
	}
	free(Y);
	Moments=Moments/G00();
	return Moments;
}


double zernike::beta()
{
	double * X = getX();
	double * Y = getY() ;
	double x_c = COF_X();
	double y_c = COF_Y();
	double A00 = 0.0;
	for(int i=0;i<M;++i)
	{
		for(int j=0;j<N;++j)
		{
			A00 += f[i][j];
		}
	}
	A00 = A00/3.14259265;
	double beta_temp = sqrt( 1/A00);
	free(X);
	free(Y);
	//cout<< beta_temp << endl;
	return beta_temp;
}


//Rho -> go to the normalized POLAR  coordinates (RHO AND THETA)
double**zernike:: RHO()
{
	double ** rho = (double**) malloc (M * sizeof(double *));
	for(int k=0;k<M;k++)
	rho[k] = (double*) malloc(N*sizeof(double));
	double * X = getX();
	double * Y = getY() ;
	double x_c = COF_X();
	double y_c = COF_Y();

	for(int i =0; i<M;i++)
	{
		for(int j=0;j<N;j++)
		{
			rho[i][j] = sqrt((X[j]-x_c)*(X[j]-x_c) + (Y[i]-y_c)*(Y[i]-y_c) );
		}
		
	}
	double sup = 0.0;
	for(int i =0; i<M;i++)
	{
		for(int j=0;j<N;j++)
		{
			if(((int)f[i][j]>0) && sup<rho[i][j])
			{
				sup=rho[i][j];
			}
		}
	}
  // cout<<sup<<endl;
	double temp_beta = beta();
	for(int i =0; i<M;i++)
	{
		for(int j=0;j<N;j++)
		{
			rho[i][j] = rho[i][j]/sup;
		
		}
	}
	free(X);
	free(Y);
	return rho;
}

double**zernike:: THETA()
{
	double ** theta = (double**) malloc (M * sizeof(double *));
	for(int k=0;k<M;k++)
		theta[k] = (double*) malloc(N*sizeof(double));
	double * X = getX();
	double * Y = getY() ;
	double x_c = COF_X();
	double y_c = COF_Y();

	for(int i =0; i<M;i++)
	{
		for(int j=0;j<N;j++)
		{
			if((X[j]-x_c)>0)
				theta[i][j] = atan((Y[i]-y_c)/(X[j]-x_c));    
			else
				theta[i][j] = 3.14259265 + atan((Y[i]-y_c)/(X[j]-x_c));  

		}
	}
	free(X);
	free(Y);
	return theta;
}

//to get factorial of integers
double zernike::factorial (int a)
{
	if (a > 1)
	{
		return (double)(a * factorial (a-1));
	}
	else
	{
		return (1.0);
	}
}

//TO CALCULATE THE CO-OFFICIENTS FOR THE ZERNIKE PLOYNOMIALS
double zernike::B(int p, int q, int k)
{
	int t = (p - abs(q))/2;
	int l = (p + abs(q))/2;
    return (  (pow( (-1.0) ,k)*(double)factorial(p-k))   /  ((double) (factorial(k)* factorial(l-k)*factorial(t-k) ))); 
	    
}

//to get the gernike moments of order p and repetition q
double* zernike::CalculateZernike(int p, int q )
{
	double * V = (double*) malloc (2* sizeof(double ));
	V[0]=0.0;V[1]=0.0;
	double** ro = RHO() ;
	double** theta = THETA();	
	int s = (p-abs(q))/2;
	for(int i =0; i<M;i++)
	{
		for(int j=0;j<N;j++)
		{
			double temp = 0.0;
			if(ro[i][j]<=1.0)
			{
				for(int k = 0 ; k<=s;k++)
				{
					temp += B(p,q,k)*pow(ro[i][j],p-2*k);
				}
				V[0] +=temp*cos(theta[i][j]*q)*f[i][j];
				V[1] -=temp*sin(theta[i][j]*q)*f[i][j];
			}
			else
			{
				V[0] += 0.0;
				V[0] += 0.0;
			}
		}
	}
	double  b = G00();
	V[0] = V[0]*(p+1)/(3.14259265*b);
	V[1] = V[1]*(p+1)/(3.14259265*b);
	free(ro);
	free(theta);
	return V;
}
std::vector< std::vector<double> > zernike::GetZernike()
{
	int p = orderofmoments;
	
	std::vector< std::vector<double> > zern_i;
	for(int i=0; i<=p; ++i)
	{
		std::vector<double> zern_j;
		for (int j=0; j<=i; ++j)
		{
			if((i-j)%2==0)
			{
				double* Z = CalculateZernike(i,j);
				if (sqrt(Z[0]*Z[0]+Z[1]*Z[1])<0.0001)
					zern_j.push_back(0);
					//cout<<"Zernike of order ("<<i<<"," <<j<<") are "<< 0 << endl;
				else
					zern_j.push_back(sqrt(Z[0]*Z[0]+Z[1]*Z[1]));
					//cout<<"Zernike of order ("<<i<<"," <<j<<") are "<<sqrt(Z[0]*Z[0]+Z[1]*Z[1])<< endl; 
				free(Z);

			}
		}
		zern_i.push_back(zern_j);
	} 
	return zern_i;
}






//
//
//
//#include <iostream>
//#include <math.h>
//#include "zernike.h"
//
//using namespace std;
///////////////////////////
///*THE CLASS ZERNIKE RETURNS THE ZERNIKE MOMENTS FOR GRAY SCALE IMAGES
//THE ZERNIKE MOMENTS ARE CALCULATED AS A LINEAR COMBINATION OF RADIAL MOMENTS
//radial moments in turn can be found as linear combinations of central geometrical moments
//Now,zernike moments are orthogonal on the unit circle(r==1) for this reason the image is mapped to the unit circle first
//then the accurate central moments */
//////////////////////////////
//
//zernike::zernike()
//{
//	M =0;
//	N =0;
//
//}
//
//zernike::~zernike()
//{
//
//}
//
///* this method is used to set the input image parameters using  Image Region Iterators
//as we need image dimensions (M and N)  for the calculations of Zernike moments and need to access the image pixel data */
//
//void zernike::SetInputImage(ImageType::Pointer inputimage)
//{
//
//	IteratorType It( inputimage, inputimage->GetRequestedRegion());
//	InputSizeType size = inputimage->GetRequestedRegion().GetSize();
//	M=size[1];
//	N=size[0];
//	ImageType::IndexType requestedIndex =inputimage->GetRequestedRegion().GetIndex();
//
//	f = (unsigned char **) malloc (M * sizeof(unsigned char *));
//	for (int k = 0; k < M; k++)
//		f[k] = (unsigned char *) malloc (N * sizeof(unsigned char));
//
//	std::cout << "M = " << M << " N = " << N << std::endl;
//
//	for (It.GoToBegin(); !It.IsAtEnd(); ++It)
//	{
//
//		ImageType::IndexType idx = It.GetIndex();			
//		f[idx[1]][idx[0]]=It.Get();
//		/*if (It.Get() >0)
//		f[idx[1]][idx[0]]=1;*/
//
//	}		
//
//}
//
//
////Now we need to map the image to the unit disk for this we are using a transformation given in the following two methods
///* this methog gives the values of X coordinates in the mapped image*/
//double* zernike::getX()
//{
//	double* X = (double *) malloc(N * sizeof(double));
//	for(int i=0;i<=N-1;++i)
//	{
//		X[i]=((double)(2*(i+1)-N-1))/(N*sqrt((double)2));
//	}
//	return X;
//}
///* this methog gives the values of Y coordinates in the mapped image*/
//double* zernike::getY()
//{
//	double* Y = (double *) malloc(M * sizeof(double));	 
//	for(int i=0;i<=M-1;++i)
//	{
//		Y[i]=(static_cast<double>(2*(i+1)-M-1))/(M*sqrt((double)2));
//		Y[i]=-1*Y[i];
//	}
//	return Y;
//}
//
//
//
///*this method returns the (0,0) order geometrical moment
//which is equal to sum(sum(f(x,y)*delX*delY);
//where delX and delY are distances between the nearest pixels in X and Y directions respectively
//*/
//double zernike::G00()
//{
//	double delY=sqrt((double)2)/M;
//	double delX=sqrt((double)2)/N;
//	double temp =0.0;
//	for(int i=0;i<M;++i)
//	{
//		for(int j=0;j<N;++j)
//		{
//			temp+=f[i][j];
//
//		}
//	}
//	//cout<<(temp*delY*delX)<<" ";
//	return (temp*delY*delX);
//}
///*this method returns the (1,0) order geometrical moment which is needed to calculate the equivalent of centre of mass for the image
//( i.e. X component of  centre of mass in image is is G10()/G00()); */
//
//double zernike::G10()
//{
//	double delY=sqrt((double)2)/M;
//	double delX=sqrt((double)2)/N;
//	double Moments=0.0;
//	double* X =getX();              
//	for(int i=0;i<M;++i)
//	{ 
//		for(int j=0;j<N;++j)
//		{
//			Moments+=((f[i][j])*X[j]);
//		}
//	}
//	free(X);
//	
//	Moments=Moments*delX*delY/G00(); 
//	return Moments;
//}
///*this method returns the (0,1) order geometrical moment which is needed to calculate the equivalent of centre of mass for the image
//( i.e. Y component of  centre of mass in image is is G01()/G00()); */
//
//
//double zernike::G01()
//{
//	double delY=sqrt((double)2)/M;
//	double delX=sqrt((double)2)/N;
//	double Moments=0.0;
//	double* Y = getY();          
//	for(int i=0;i<M;++i)
//	{ 
//		for(int j=0;j<N;++j)
//		{
//			Moments+=((f[i][j])*Y[i]);
//		}
//	}
//	free(Y);
//
//	Moments=Moments*delX*delY/G00(); 
//	return Moments;
//}
///* Ip is an array used to generate the integral of (x^p.delx) which is needed in the accurate calculation of geometrical moments */
//
//double* zernike::generateIp(int p)
//{   
//	double* X  =  getX();
//	double* Ip = (double *) malloc(N * sizeof(double));
//	double U1  = 0.0,U2=0.0;
//	double m0  = (G10());
//	double del = sqrt((double)2)/(N*2.0);
//
//	for(int i=0;i<N;++i)
//	{
//		U1 = X[i] -m0+del;
//		U2 = X[i] -m0-del;
//		Ip[i] = (1.0/(double)(p+1))*( pow((U1),(p+1)) - pow(U2,(p+1)));
//	}
//	free(X);
//	
//	return Ip;
//}
//
///*Iq is an array used to generate the integral of (y^p.dely) which is needed in the accurate calculation of geometrical moments */
//double* zernike::generateIq(int q)
//{  
//	double* Y  = getY();
//	double* Iq = (double *) malloc(M * sizeof(double));
//	double V1  = 0.0,V2=0.0;
//	double m   = G01();
//	double del = sqrt((double)2)/(M*2.0);
//	for(int i=0;i<M;++i)
//	{
//		V1 = Y[i] -m+del;
//		V2 = Y[i] - m-del;
//		Iq[i] = (1.0/(double)(q+1))*( pow((V1),(q+1)) - pow(V2,(q+1)));
//	}
//	free(Y);
//	
//	return Iq;
//}
//
//
////to  generate the geometrical moments 
//double zernike::generateGmoments(int p,int q)	
//{ 
//	double Moments;
//	Moments=0; 
//	double* Ip = generateIp(p);
//	double* Iq = generateIq(q);
//	
//    for(int i=0;i<M;++i)
//		{
//			for(int j=0;j<N;++j)
//			{
//				Moments=Moments+(Ip[j]*Iq[i]*f[i][j]);
//
//			}
//		}
//	/*for(int i=0;i<M;++i)
//	{
//		Moments=Moments+Iq[i]*Yq[i];
//	}
//	for(int i=0;i<M;++i)
//	{
//		for(int j=0;j<N;++j)
//		{
//			Moments=Moments+(Ip[j]*Iq[i]*f[i][j]);
//			
//		}
//	}*/
//	/*double* Yq =(double *) malloc(M* sizeof(double));
//	for(int i=0;i<M;++i)
//	{
//		Yq[i]=0.0;
//	}
//
//	for(int i=0;i<M;++i)
//	{
//		for(int j=0;j<N;++j)
//		{
//			Yq[i] = Yq[i]+(Ip[j]*f[i][j]);
//			
//		}
//	}
//
//
//
//	if(flag)
//	{
//		for(int i=0;i<M;++i)
//		{
//			for(int j=0;j<N;++j)
//			{
//				Moments=Moments+(Ip[j]*Iq[i]*f[i][j]);
//
//			}
//		}
//	}
//	else
//	{
//		for(int i=0;i<M;++i)
//		{
//			for(int j=0;j<N;++j)
//			{
//				if (f[i][j]>0)
//				{
//					f[i][j]=255;
//					Moments=Moments+(Ip[j]*Iq[i]*f[i][j]);
//				}
//			}
//		}
//	}*/
//
//	free(Ip);
//	free(Iq);
//	//free(Yq);
//
//	return Moments; 
//}
//
////to get factorial of integers
//double zernike::factorial (int a)
//{
//	if (a > 1)
//	{
//		return (double)(a * factorial (a-1));
//	}
//	else
//	{
//		return (1.0);
//	}
//}
//
//
////to get combination of two integers
//double zernike::combination(int m,int n)
//{
//	if(n==0)
//	{
//		return 1;
//	}
//	else if(m == n)
//	{
//		return 1;
//	}
//	else if (m>n)
//	{
//		return ((double)((static_cast<double> (m))/(static_cast<double>(m-n)))*combination(m-1,n));
//	}
//	else
//	{
//		cout << "ERROR: Combination failed: m = " << m << " n = " << n << endl;
//		exit(0);
//	}
//}
//
//
//double zernike::cosine1(int m) 
//{
//	double x=cos(m*3.14159265/2.0);
//	if((x<0.99) && (x >(-0.99)) )
//		return 0.0;
//	else if((x<=(-0.99)))
//		return -(1.0);
//	else
//		return 1.0;
//}
//
//double zernike::sine1(int m)
//{  
//	double y=sin(m*3.14159265/2.0);
//	if(((y<0.99 )&& (y>(-0.99))))
//		return 0.0;
//	else if((y<=(-0.99)))
//		return -(1.0);
//	else
//		return 1.0;
//}
////to calculate the radial moments .to calculate radial moments the geometrical moments are normalized by pow(generateGmoments(0,0),p/2+1) 
//double* zernike::Radialmoments(int p,int q)
//{  
//	double temp=0.0;
//	double G=generateGmoments(0,0);
//	double scalefactor=(pow(sqrt(G),p))*G;
//	double RadReal=0,RadIm=0;                         // for real and imaginary part of RadialMoments
//	double* R = (double *) malloc(2 * sizeof(double));
//
//	int s;
//	s=(p-q)/2;
//	for(int j=0;j<=s;++j)
//	{
//		for(int m=0;m<=q;++m)
//		{
//			temp     = combination(s,j)* combination(q,m)*generateGmoments((p-q-2*j+m),(2*j-m+q))/scalefactor;
//			
//			RadReal += (cosine1(m)*temp);
//			RadIm   += ((-1.0)*sine1(m)  * temp);
//		}
//	}
//	R[0]  =  RadReal; 
//	R[1]  =  RadIm;
//	return R;
//}
//	
//
///* as mentioned above the zernike moments are calculated as a linear combination of the radial moments.Bpqk are the multiplying coefficients */
//
//double zernike::B(int p, int q, int k)
//{
//	//return (  (pow( (-1.0) ,(p-k)/2)*(double)factorial((p+k)/2))   /  ((double) (factorial((p-k)/2)* (double)factorial((k+q)/2)*(double)factorial((k-q)/2)) )); 
//	
//	if(p==q && q==k)
//		return 1.0;
//	else if(q!=p && k==p)
//		return (((p+q+2)/(double)(p-q))*(B(p,q+2,p)));
//	else
//		return ((double)(k+2+q)*(double)(k+2-q)/((double)(p+k+2)*(double)(p-k))*B(p,q+2,p));
//
//}
//
//// Calculation of Zernike Moments Using Radial Moments.
//
//double* zernike::CalculateZernike(int p, int q)
//{
//	
//	int k=0;
//	double pi = 3.14259265;
//	double* Z = (double *) malloc(2 * sizeof(double));
//
//	for (int i = 0; i < 2; i++)
//		Z[i] = 0;
//
//	if (p>=abs(q) && (p-q)%2==0)		
//	{
//		for(k=abs(q); k<=p ; k+=2)
//		{
//			double *R = Radialmoments(k,abs(q));
//			double b = B(p,abs(q),k);
//			Z[0] =Z[0]+( b * R[0]);
//			Z[1] =Z[1]+( b * R[1]);
//			free(R);
//			
//		}
//		Z[0] =Z[0]*(double)(p+1)/pi;
//		Z[1] =Z[1]*(double)(p+1)/pi;
//
//		if(q<=0)
//		{
//			Z[0] = Z[0];
//			Z[1] = -1*Z[1];
//			return Z;			
//		}		
//		else
//			return Z;
//	}
//
//	
//	else
//	{
//		cout << "Error: p is less than q or p-q is odd: p = " << p << " q = " << q << endl; 
//		exit(0);
//	}
//}
//void zernike::GetZernike(ImageType::Pointer inputimage,int p)
//{
//	SetInputImage(inputimage);
//	for(int i=0; i<=p; ++i)
//	{
//		for (int j=0; j<=i; ++j)
//		{
//			if((i-j)%2==0)
//			{
//
//				double* Z = CalculateZernike(i,j);
//				if (sqrt(Z[0]*Z[0]+Z[1]*Z[1])<0.00001)
//					cout<<"Zernike of order ("<<i<<"," <<j<<") are "<< 0 << endl;
//				else
//					cout<<"Zernike of order ("<<i<<"," <<j<<") are "<<sqrt(Z[0]*Z[0]+Z[1]*Z[1])<< endl; 
//				free(Z);
//			}
//				
//			
//		} 
//	}  
//}

} //end namespace ftk
