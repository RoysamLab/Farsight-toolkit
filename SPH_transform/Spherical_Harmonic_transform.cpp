#include "Spherical_Harmonic_transform.h"

#ifdef _OPENMP
#include "omp.h"
#endif

SPH_Trasform::SPH_Trasform()
{

}

SPH_Trasform::~SPH_Trasform()
{
}

void SPH_Trasform::Initialize( vnl_vector<double> x, vnl_vector<double> y, vnl_vector<double> z)
{
	if( x.size() != y.size() || y.size() != z.size() || z.size() != x.size() )
	{ 
		std::cout<<"invalida data, unconsistent size of vector x, y, z !"<<std::endl;
		return;
	}
	this->x.clear();
	this->y.clear();
	this->z.clear();
	this->R.clear();
	this->lat.clear();
	this->lon.clear();
	this->theta.clear();
	this->phi.clear();
	this->fphth.clear();
	this->lmcosi.clear();
	this->coslat.clear();
	this->plm.clear();
	this->C.clear();

	this->x = x;
	this->y = y;
	this->z = z;
	this->numpoints = x.size();
}

void SPH_Trasform::Initialize( std::vector<std::vector<double> > data, int L)
{
	if( data.empty () || data[0].empty () )
	{ 
		std::cout<<"invalida data, unconsistent size of vector x, y, z !"<<std::endl;
		return;
	}
	this->x.clear();
	this->y.clear();
	this->z.clear();
	this->R.clear();
	this->lat.clear();
	this->lon.clear();
	this->theta.clear();
	this->phi.clear();
	this->fphth.clear();
	this->lmcosi.clear();
	this->coslat.clear();
	this->plm.clear();
	this->C.clear();

	this->numpoints = data.size();
	this->x.set_size(this->numpoints);
	this->y.set_size(this->numpoints);
	this->z.set_size(this->numpoints);
	for(unsigned int i = 0; i < this->numpoints; i++)
	{
		this->x(i) = data[i][0];
		this->y(i) = data[i][1];
		this->z(i) = data[i][2];
	}
	this->L = L;
}

void SPH_Trasform::xyz2rll()
{
	for(unsigned int i = 0; i < this->numpoints; i++)
	{
		if ( abs ( this->x(i) - this->y(i) ) < EPS && abs ( this->x(i) ) < EPS )
		{
			this->x(i) = 0.01;
			this->y(i) = 0.01;
			std::cout<<"there is coordinate change !"<<std::endl;
		}
	}

	this->R.set_size(this->numpoints);
	this->theta.set_size(this->numpoints);
	this->phi.set_size(this->numpoints);

	for(unsigned int i = 0; i < this->numpoints; i++)
	{
		this->R(i) = sqrt( this->x(i) * this->x(i) + this->y(i) * this->y(i) +this->z(i) * this->z(i));
		this->theta(i) = acos( this->z(i) / this->R(i) );
		this->phi(i) = acos( this->x(i) / sqrt( this->x(i) * this->x(i) + this->y(i) * this->y(i)) );
		if(this->y(i) < 0)
		{
			this->phi(i) += PI;
		}
	}
}

void SPH_Trasform::DataFilter()
{
	this->fphth.set_size(this->numpoints);
	this->lat.set_size(this->numpoints);
	this->lon.set_size(this->numpoints);

	//Sort ascending acording to theta
	std::vector<int > index;
	for (unsigned int i = 0; i < this->theta.size () ; i++)
		index.push_back (i);
	SortwithIndex(this->theta, 0, this->theta.size() - 1, index);
	for(int i = 0; i < index.size(); i++)
	{
		this->fphth(i) = this->R(index[i]);
		this->lon(i) = this->phi(index[i]);
		this->lat(i) = this->theta(i);
	}

	//////////////////////////////////for debuging/////////////////////////////////////
	//const char* filename2 = "justsorting.txt";
	//FILE *fp2 = fopen(filename2,"w");
	//fprintf(fp2,"%s", "fthph");
	//fprintf(fp2,"\t");
	//fprintf(fp2,"%s", "lon");
	//fprintf(fp2,"\t");
	//fprintf(fp2,"%s", "lat");
	//fprintf(fp2,"\n");

	//for(int i = 0; i < this->numpoints; i++)
	//{
	//	fprintf( fp2,"%f", this->fphth(i) );
	//	fprintf(fp2,"\t");
	//	fprintf( fp2,"%f", this->lon(i) );
	//	fprintf(fp2,"\t");
	//	fprintf(fp2,"%f", this->lat(i) );
	//	fprintf(fp2,"\n");
	//}
	//fclose(fp2);
	//////////////////////////////////for debuging/////////////////////////////////////

	//sorting phi within same theta in order
	unsigned int i = 0;
	while (i < this->numpoints - 1)
	{
		unsigned int j = 1;
		while (i + j < this->numpoints)
		{
			if ( abs ( this->lat(i) - this->lat(i + j) )< EPS ) 
				j = j + 1;
			else
				break;
		}
		if ( j > 1)
		{
			std::vector<int > tempindex;
			for(unsigned int k = 0; k < this->lon.size(); k++ )
				tempindex.push_back (k);
			SortwithIndex( this->lon, i, i + j - 1, tempindex ); 

			std::vector<double > tempfthph;
			for(unsigned int p = 0;  p < j;  p++)
				tempfthph.push_back ( this->fphth(p + i) );
			for(unsigned int p = 0; p < j; p++)
				this->fphth( p + i ) = tempfthph[ tempindex[p + i] - i ];
		}
		i = i + j;
	}
	//////////////////////////////////for debuging/////////////////////////////////////
	//const char* filename1 = "againsorting.txt";
	//FILE *fp1 = fopen(filename1,"w");
	//fprintf(fp1,"%s", "fthph");
	//fprintf(fp1,"\t");
	//fprintf(fp1,"%s", "lon");
	//fprintf(fp1,"\t");
	//fprintf(fp1,"%s", "lat");
	//fprintf(fp1,"\n");

	//for(int i = 0; i < this->numpoints; i++)
	//{
	//	fprintf( fp1,"%f", this->fphth(i) );
	//	fprintf(fp1,"\t");
	//	fprintf( fp1,"%f", this->lon(i) );
	//	fprintf(fp1,"\t");
	//	fprintf(fp1,"%f", this->lat(i) );
	//	fprintf(fp1,"\n");
	//}
	//fclose(fp1);
	//////////////////////////////////for debuging/////////////////////////////////////

	//check and delete the same coordinate
	i = 0;
	int k = 0;
	while (i < this->numpoints - 1)
	{
		unsigned int j = 1;
		while (i + j < this->numpoints)
		{
			if ( abs( this->lat(i) - this->lat(i + j) ) < EPS )
			{
				if ( abs( this->lon(i) - this->lon(i + j) ) < EPS)
				{
					if ( abs( this->fphth(i) - this->fphth(i + j) ) < EPS )
					{
						j = j + 1;
						std::cout<<"there is same coordinate, deleted while converting !"<<std::endl;
					}
					else
						break;
				}
				else
					break;
			}
			else
				break;
		}
		this->fphth(k) = this-> fphth(i);
		this->lat(k) = this-> lat(i);
		this->lon(k) = this-> lon(i);
		k = k + 1;
		i = i + j;    
	}
	if( i == this->numpoints - 1)
	{
		this->fphth(k) = this-> fphth(i);
		this->lat(k) = this-> lat(i);
		this->lon(k) = this-> lon(i);
		k = k + 1;
	}
	this->numcoord = k;

	//////////////////////////////////for debuging/////////////////////////////////////
	//const char* filename = "justdelete.txt";
	//FILE *fp = fopen(filename,"w");
	//fprintf(fp,"%s", "fthph");
	//fprintf(fp,"\t");
	//fprintf(fp,"%s", "lon");
	//fprintf(fp,"\t");
	//fprintf(fp,"%s", "lat");
	//fprintf(fp,"\n");

	//for(int i = 0; i < this->numcoord; i++)
	//{
	//	fprintf( fp,"%f", this->fphth(i) );
	//	fprintf(fp,"\t");
	//	fprintf( fp,"%f", this->lon(i) );
	//	fprintf(fp,"\t");
	//	fprintf(fp,"%f", this->lat(i) );
	//	fprintf(fp,"\n");
	//}
	//fclose(fp);
	//////////////////////////////////for debuging/////////////////////////////////////
}

void SPH_Trasform::SortwithIndex( vnl_vector<double> &elearray, int start, int end, std::vector<int > &index)
{
	int pointer1 = start;
	int pointer2 = end;
	double pivotvalue = elearray(start);
	int indexstart = index[start];

	while( pointer1 < pointer2 )
	{
		while( (elearray(pointer2) >= pivotvalue) && (pointer1 < pointer2) )
			pointer2--;	
		if (pointer1 != pointer2)
		{
			index[pointer1] = index[pointer2];
			elearray(pointer1) = elearray(pointer2);
			pointer1++;
		}
		while ( (elearray(pointer1) <= pivotvalue) && (pointer1 < pointer2) )
			pointer1++;

		if (pointer1 != pointer2)
		{
			index[pointer2] = index[pointer1];
			elearray(pointer2) = elearray(pointer1);
			pointer2--;
		}
	}
	elearray(pointer1) = pivotvalue;
	index[pointer1] = indexstart;
	int pointer = pointer1;

	int i = 0;
	if(start < pointer)
	{
		SortwithIndex(elearray, start, pointer-1, index);
	}
	if(end > pointer)
	{
		SortwithIndex(elearray, pointer+1, end, index);
	}
}

void SPH_Trasform::costheta()
{
	this->coslat.set_size( this->numcoord );
	for (unsigned int i = 0; i < this->numcoord; i++)
		this->coslat(i) = cos ( this->lat(i) );
	//////////////////////////////////for debuging/////////////////////////////////////
	//const char* filename = "cosine.txt";
	//FILE *fp = fopen(filename,"w");

	//for(int i = 0; i < this->numcoord; i++)
	//{
	//	fprintf( fp,"%f", this->coslat(i) );
	//	fprintf(fp,"\n");
	//}
	//fclose(fp);
	//////////////////////////////////for debuging/////////////////////////////////////
}

void SPH_Trasform::Legendre_P()
{
	this->plm.set_size( this->numcoord, this-> size_coefficient * 2);
	for(unsigned int k = 0; k < this->numcoord; k++ )
	{ 
		int counter = 0;
		for(int l = 0; l <= this->L; l++)
		{
			int innercounter = 0;
			for(int m = 0; m <= l; m++)
			{
				double fractial = Fractial(l, m);
				double factor = sqrt((double)(2*l+1));
				double l_p = legendre(l, m, this->coslat(k) ) ;
				this->plm(k, counter + innercounter) =  l_p * fractial * factor;
				if(m > 0)
					this->plm(k, counter + innercounter) *= SQROOT;
				if((m % 2) != 0)
					this->plm(k, counter + innercounter) *= (-1);
				innercounter++;
			}
			counter += l + 1;
		}
	}
	//////////////////////////////////for debuging/////////////////////////////////////
	//const char* filename = "plm.txt";
	//FILE *fp = fopen(filename,"w");

	//for(int i = 0; i < this->numcoord; i++)
	//{
	//	for(int j = 0 ;j < this->size_coefficient; j++)
	//	{
	//		fprintf( fp,"%f", this->plm(i,j) );
	//		fprintf(fp,"\t");
	//	}
	//	fprintf(fp,"\n");
	//}
	//fclose(fp);
	//////////////////////////////////for debuging/////////////////////////////////////
}

void SPH_Trasform::GetCoefficientSize()
{
	int counter = 0;
	for(int i = 0; i <= L; i++)
	{
		counter += i + 1;
	}
	this->size_coefficient = counter;
}

double SPH_Trasform::legendre( int l, int m, double x )
{
	double le_p;
	le_p = boost::math::legendre_p(l, m, x);
	return le_p;
}

void SPH_Trasform::CombiningLon()
{
	//creat m vector for folling steps
	std::vector<int > m;
	for(int i = 0; i <= L; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			m.push_back(j);
		}
	}

	//calculate cosine and sine binding of plm
	vnl_matrix<double> cosmat;
	cosmat.set_size( this->numcoord, this->size_coefficient );
	vnl_matrix<double> sinmat;
	sinmat.set_size( this->numcoord, this->size_coefficient );
	for(unsigned int i = 0; i < this->numcoord; i++)
	{
		for(unsigned int j = 0; j < this->size_coefficient;  j++)
		{
			double temp =  this->lon(i) * m[j] ;
			cosmat(i,j) = cos ( temp );
			sinmat(i,j) = sin ( temp ) ;
		}
	}

	//cosine and sinine binding of plm
	for(unsigned int i = 0; i < this->numcoord; i++)
	{
		for(unsigned int j = 0; j < this->size_coefficient;  j++)
		{
			this->plm(i,j + size_coefficient) = sinmat(i,j) * this->plm(i,j);
		}
	}
	for(unsigned int i = 0; i < this->numcoord; i++)
	{
		for(unsigned int j = 0; j < this->size_coefficient;  j++)
		{
			this->plm(i,j) = cosmat(i,j) * this->plm(i,j);
		}
	}

	//////////////////////////////////for debuging/////////////////////////////////////
	const char* filename = "plmcombining.txt";
	FILE *fp = fopen(filename,"w");

	for(unsigned int i = 0; i < this->numcoord; i++)
	{
		for(int j = 0 ;j < this->size_coefficient * 2; j++)
		{
			fprintf( fp,"%f", this->plm(i,j) );
			fprintf(fp,"\t");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	/////////////////////////////////////////////////////////
	//const char* filename1 = "plmcos.txt";
	//FILE *fp1 = fopen(filename1,"w");

	//for(unsigned int i = 0; i < this->numcoord; i++)
	//{
	//	for(int j = 0 ;j < this->size_coefficient; j++)
	//	{
	//		fprintf( fp1,"%f", cosmat(i,j) );
	//		fprintf(fp1,"\t");
	//	}
	//	fprintf(fp1,"\n");
	//}
	//fclose(fp1);
	//////////////////////////////////for debuging/////////////////////////////////////
}

void SPH_Trasform::DataFitting()
{
	std::cout<<"Calculating SVD !"<<std::endl ;
	vnl_svd< double > svd(this->plm);
	vnl_matrix< double > U = svd.U();
	vnl_diag_matrix< double > W = svd.W();
	vnl_matrix< double > V = svd.V();
	std::cout<<"Done Calculating SVD !"<<std::endl ;
	//////////////////////////////////for debuging/////////////////////////////////////
	const char* filename = "U.txt";
	FILE *fp = fopen(filename,"w");

	for(unsigned int i = 0; i < U.rows() ; i++)
	{
		for(int j = 0 ;j <U.columns(); j++)
		{
			fprintf( fp,"%f", U(i,j) );
			fprintf(fp,"\t");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	const char* filename1 = "V.txt";
	FILE *fp1 = fopen(filename1,"w");

	for(unsigned int i = 0; i < V.rows() ; i++)
	{
		for(int j = 0 ;j <V.columns(); j++)
		{
			fprintf( fp1,"%f", V(i,j) );
			fprintf(fp1,"\t");
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);

	const char* filename2 = "w.txt";
	FILE *fp2 = fopen(filename2,"w");

	for(unsigned int i = 0; i < W.rows() ; i++)
	{
			fprintf( fp2,"%f", W(i,i) );
			fprintf(fp2,"\n");
	}
	fclose(fp2);
	////////////////////////////////////////////////////for debuging////////////////////////////////////////////////
	for(unsigned int i = 0; i < W.rows (); i++)
	{
		if (abs(W(i,i)) < EPS)
			W(i,i) = 0;
		else
			W(i,i) = 1 / W(i,i);
	}

	vnl_matrix< double > firstproduct;
	firstproduct.set_size(V.rows (), W.columns () );
	vnl_matrix< double > sencondproduct;
	sencondproduct.set_size(V.rows (), U.rows () );
	if( this->numcoord > this->size_coefficient * 2 )
	{
		firstproduct = V * W;
		sencondproduct = firstproduct * U.transpose () ;
	}
	else
	{
		std::cout<<"there is less equa than unknows !"<<std::endl;
	}

	vnl_vector<double> tempfphth;
	tempfphth.set_size(this->numcoord);
	for(unsigned int i = 0; i < this->numcoord; i++)
		tempfphth(i) = this->fphth(i);
	this->C.set_size(this->size_coefficient * 2);
	this->C = tempfphth.pre_multiply(sencondproduct);
}

void SPH_Trasform::lmCalculating()
{
	this->lmcosi.set_size( this->size_coefficient, 4 );
	int counter = 0;
	for(int l = 0; l <= L; l++)
	{
		for(int j = 0; j <= l; j++)
		{
			this->lmcosi( counter,0 ) = (double)l;
			this->lmcosi( counter,1 ) = (double)j;
			counter++;
		}
	}

	for(unsigned int i = 0; i < this->size_coefficient; i++ )
	{
		this->lmcosi( i,2 ) = this->C(i);
		this->lmcosi( i,3 ) = this->C(i + this->size_coefficient );
	}
	//////////////////////////////////for debuging/////////////////////////////////////
	const char* filename = "whatthehell.txt";
	FILE *fp = fopen(filename,"w");

	for(unsigned int i = 0; i < this->size_coefficient; i++)
	{
			fprintf( fp,"%f", this->lmcosi(i,0) );
			fprintf(fp,"\t");
			fprintf( fp,"%f", this->lmcosi(i,1) );
			fprintf(fp,"\t");
			fprintf( fp,"%f", this->lmcosi(i,2) );
			fprintf(fp,"\t");
			fprintf( fp,"%f", this->lmcosi(i,3) );
			fprintf(fp,"\n");
	}
	fclose(fp);
	/////////////////////////////////////////////////////////
}

void SPH_Trasform::Transform()
{
	this->xyz2rll ();
	this->DataFilter ();
	this->costheta ();
	this->GetCoefficientSize ();
	this->Legendre_P ();
	this->CombiningLon ();
	this->DataFitting();
	this->lmCalculating ();
}

void SPH_Trasform::WriteFile()
{
	const char* filename = "lmcosi.txt";
	FILE *fp = fopen(filename,"w");
	// head names
	fprintf(fp,"%s", "degree");
	fprintf(fp,"\t");
	fprintf(fp,"%s", "order");
	fprintf(fp,"\t");
	fprintf(fp,"%s", "cosine");
	fprintf(fp,"\t");
	fprintf(fp,"%s", "sine");
	fprintf(fp,"\n");

	for(unsigned int i = 0; i < this->size_coefficient; i++)
	{
		fprintf( fp,"%d", this->lmcosi(i,0) );
		fprintf(fp,"\t");
		fprintf( fp,"%d", this->lmcosi(i,1) );
		fprintf(fp,"\t");
		fprintf(fp,"%f", this->lmcosi(i,2) );
		fprintf(fp,"\t");
		fprintf(fp,"%f", this->lmcosi(i,3) );
		fprintf(fp,"\n");
	}
	fclose(fp);
}

double SPH_Trasform::Fractial(int l, int m)
{
	size_t fractial = 1;
	int num = l - m;
	int den = l + m;
	if(den - num > 15)
		return 0.0;
	while(den > num)
	{
		fractial *= den; //fractial = fractial * den
		den--;
	}

	double fractial_reciprocal = (double)(1 / (double)fractial);

	//assert(isValid(fractial_reciprocal)); //We did not get a valid number if this fails
		
	return sqrt (fractial_reciprocal);
}