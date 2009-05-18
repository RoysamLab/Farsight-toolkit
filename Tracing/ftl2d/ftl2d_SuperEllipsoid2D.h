#ifndef SUPERELLIPSOID2D_H
#define SUPERELLIPSOID2D_H

#include <itkImage.h>
#include <itkSmartPointer.h>
#include <itkImageRegionIterator.h> 
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageLinearConstIteratorWithIndex.h>

#include <vnl\vnl_math.h>
#include <vnl\vnl_matrix_fixed.h>
#include <vnl\vnl_vector_fixed.h>
#include <string>

typedef itk::Image< float, 2 >   ImageType;
typedef itk::Point< double, 2 >   PointType;

typedef vnl_vector_fixed <double,2> Vect2;
typedef vnl_matrix_fixed <double,2,2> Mat22;

#define NN 29


class SuperEllipsoid2D	{

public:	
	SuperEllipsoid2D(double);
	~SuperEllipsoid2D();
	
	void generateConvexHull(double&, double&, double&, Vect2&, double&);
	bool Interpolation(ImageType::Pointer);
	void UpdateLikelihood(float&, float&);
	void gradient_superquad2d(double&, double&, double&, Vect2& , double&);
	
	void UpdateOrientation(double&, unsigned int);
	void UpdatePosition(Vect2& , unsigned char , double&);
	void UpdateScale(double&, double&, double&);
	void UpdateShapeParameter(double&, double&, double&, double&);
	void UpdateAxis(double&, double&, double&, Vect2&);
	void PrintSelf();
	

private:		
	
	double step;
	
	vnl_matrix_fixed <double,NN+1,2> s;	//oriented SQ
	vnl_matrix_fixed <double,NN+1,2> S;	//scaled SQ
	vnl_matrix_fixed <double,NN+1,2> U;	//unscaled SQ
	vnl_matrix_fixed <double,NN+1,2> N;	//unscaled SQ
	vnl_vector_fixed <double,NN+1> L;		//Likelihoods
	vnl_vector_fixed <double,NN+1> A;		//area (length)
	vnl_vector_fixed <double,NN+1> F;		//image values
	double MaxAspectRatio;
	
};

#endif
