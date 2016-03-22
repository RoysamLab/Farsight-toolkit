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

void SPH_Trasform::Initialize(  std::vector<std::vector<double> > data, int L)
{
	if( data.empty() || data[0].size() != 3)
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

	for(int i = 0; i < data.size(); i++)
	{
		this->x(i) = data[i][0];
		this->y(i) = data[i][1];
		this->z(i) = data[i][z];
	}
	this->numpoints = data.size();
	this->L = L;
}

void SPH_Trasform::xyz2rll()
{
	this->R.set_size(this->numpoints);
	this->theta.set_size(this->numpoints);
	this->phi.set_size(this->numpoints);

	for(int i = 0; i < this->numpoints; i++)
	{
		this->R(i) = sqrt( this->x(i) * this->x(i) + this->y(i) * this->y(i) +this->z(i) * this->z(i));
		this->theta(i) = acos( this->z(i) / this->R(i) );
		this->phi(i) = acos( this->x(i) / sqrt( this->x(i) * this->x(i) + this->y(i) * this->y(i)) );
	}
}

void SPH_Trasform::DataFilter()
{
	this->fphth.set_size(this->numpoints);
	this->lat.set_size(this->numpoints);
	this->lon.set_size(this->numpoints);

	//Sort ascending acording to theta
	std::vector<int > index = SortwithIndex(this->theta, 0, this->theta.size() - 1);
	for(int i = 0; i < index.size(); i++)
	{
		this->fphth(i) = this->R(index[i]);
		this->lon(i) = this->phi(index[i]);
		this->lat(i) = this->theta(i);
	}

	//check and delete the same coordinate
	int i = 0;
	int k = 0;
	while (i < this->numpoints - 1)
	{
		int j = 1;
		while (i + j < this->numpoints)
		{
			if ( this->lat(i) == this->lat(i + j) )
			{
				if ( this->lon(i) == this->lon(i + j) )
				{
					if ( this->fphth(i) == this->fphth(i + j) )
					{
						j = j + 1;
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
	if( i = this->numpoints - 1)
	{
		this->fphth(k) = this-> fphth(i);
		this->lat(k) = this-> lat(i);
		this->lon(k) = this-> lon(i);
		k = k + 1;
	}
	this->numcoord = k;

	//rearange phi within same theta in order
	int i = 0;
	while (i < this->numcoord - 1)
	{
		int j = 1;
		while (i + j < this->numcoord)
		{
			if ( this->lat(i) == this->lat(i + j) ) 
				j = j + 1;
			else
				break;
		}
		if ( j > 1)
		{
			std::vector<int > tempindex = SortwithIndex( this->lon, i, i + j ); 
			std::vector<double > tempfthph;
			for( int p = 0; p < tempindex.size(); p++ )
				tempfthph.push_back( this->fphth( i + p ) ); 
			for( int p = 0; p < tempindex.size(); p++ )
				this->fphth( i + p ) = tempfthph[p];
		}
		i = i + j;
	}
}

std::vector<int > SPH_Trasform::SortwithIndex( vnl_vector<double> elearray, int start, int end)
{
	int pointer1 = start;
	int pointer2 = end;
	double pivotvalue = elearray(start);

	std::vector<int > index;
	for(int  i = 0; i < (end - start + 1); i++)
		index.push_back(start + i);

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
	index[pointer1] = index[start];
	int pointer = pointer1;

	int i = 0;
	if(start < pointer)
	{
		std::vector<int > index1 = SortwithIndex(elearray, start, pointer-1);
		for(; i < index1.size(); i++)
			index[i] = index1[i];
	}
	if(end > pointer)
	{
		std::vector<int > index2 = SortwithIndex(elearray, pointer+1, end);
		for(; i < index2.size(); i++)
			index[i] = index2[i];
	}
	return index;
}

void SPH_Trasform::costheta()
{
	this->coslat.set_size( this->numcoord );
	for (int i = 0; i < this->numcoord; i++)
		this->coslat(i) = cos ( this->lat(i)*PI/180 );
}

void SPH_Trasform::Legendre_P()
{
	this->plm.set_size( this->numcoord, this-> size_coefficient * 2);
	for( int k = 0; k < this->numcoord; k++ )
	{ 
		int counter = 0;
		for(int l = 0; l <= this->L; l++)
		{
			int innercounter = 0;
			for(int m = 0; m <= l; m++)
			{
				this->plm(k, counter + innercounter) = legendre(l, m, this->coslat(i) ) ;
				innercounter++;
			}
			counter += l + 1;
		}
	}
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
	le_p = boost::math::legendre_p(m, m, x);
	return le_p;
}

void SPH_Trasform::CombiningLon()
{
	//creat m vector for folling steps
	std::vector<int > m;
	m.resize(this->size_coefficient);
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
	for(int i = 0; i < this->numcoord; i++)
	{
		for(int j = 0; j < this->size_coefficien;  j++)
		{
			cosmat(i,j) = cos ( this->lon(i) * m(j) * PI / 180 );
			sinmat(i,j) = sin ( this->lon(i) * m(j) * PI / 180) * this->plm(i,j);
		}
	}

	//cosine and sinine binding of plm
	for(int i = 0; i < this->numcoord; i++)
	{
		for(int j = 0; j < this->size_coefficien;  j++)
		{
			this->plm(i,j) = cosmat(i,j) * this->plm(i,j);
		}
		for(; j < this->size_coefficien * 2 - 1;  j++)
		{
			this->plm(i,j) = sinmat(i,j);
		}
	}
}

void SPH_Trasform::DataFitting()
{
	vnl_svd< double > svd(this->plm);
	vnl_matrix< double > U = svd.U();
	vnl_diag_matrix< double > W = svd.W();
	vnl_matrix< double > V = svd.V();
	
	for(int i = 0; i < W.rows (); i++)
	{
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

	this->C.set_size(this->size_coefficient * 2);
	this->C = this->fphth.pre_multiply(sencondproduct);
}

void SPH_Trasform::lmCalculating()
{
	this->lmcosi.set_size( this->size_coefficient, 4 );
	int counter = 0;
	for(int i = 0; i <= L; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			this->lmcosi( counter,0 ) = counter;
			this->lmcosi( counter,1 ) = j;
			counter++;
		}
	}

	for( int i = 0; i < this->size_coefficient; i++ )
	{
		this->lmcosi( i,2 ) = this->C(i);
		this->lmcosi( i,3 ) = this->C(i + this->size_coefficient );
	}
}

void SPH_Trasform::Transform()
{
	this->xyz2rll();
	this->DataFilter();
	this->costheta();
	this->GetCoefficientSize();
	this->Legendre_P();
	this->CombiningLon();
	this->DataFitting();
	this->lmCalculating();
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

	//content
	for(int i = 0; i < this->size_coefficient; i++)
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