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

#ifndef SEGMENT2D_H
#define SEGMENT2D_H

#include <itkImage.h>
#include <itkSmartPointer.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageLinearConstIteratorWithIndex.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>


#include <string>

typedef itk::Image< float, 2 >   ImageType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
typedef vnl_vector_fixed <double,2> Vect2;
typedef vnl_matrix_fixed <double,2,2> Mat22;

#include "ftl2d_SeedPoint2D.h"
#include "ftl2d_TraceConfig.h"
#include "ftl2d_SuperEllipsoid2D.h"

class Segment2D	{
public:
	Segment2D(TraceConfig *);
	Segment2D(Segment2D *);
	~Segment2D();
	int getID() {return SegmentID;}
	void Initialize(ImageType::Pointer);
	void advance(Segment2D *);
	void CopySeedInformation(SeedPoint2D* seed);
	bool IsSeedSegmentValid();
	bool IsSegmentValid();
	void fitSuperEllipsoid(ImageType::Pointer);
	bool UpdateIntensity(ImageType::Pointer);
	void EstimatorMAD(std::vector<float>&, std::vector<float>&);
	void EstimatorMAD2(std::vector<float>&, std::vector<float>&);
	double ComputeMedian(std::vector<float>&);
	void flipDirection(unsigned int);
	void SegmentValidCheck();
	void PrintSelf(unsigned char);
	double getLength();
	
	inline float getForegdIntensity() {return(this->ForegdIntensity);}
	inline float getBackgdIntensity() {return(this->BackgdIntensity);}
	inline double getq() {return (this->q);}
	inline double geta1() {return (this->a1);}
	inline double geta2() {return (this->a2);}
	inline double gete1() {return (this->e1);}
	inline double getMAD() {return (this->MADEstimate);}
	inline double getLikelihood() {return (this->Likelihood);}
	inline Vect2 getmu() {return (this->mu);}
	inline unsigned int getEndType(){return(this->EndType);};
	inline Vect2 getR2(){return (this->R2);}
	inline double getPredictedStep() {return(this->PredictedStep);}
	inline void setEndType(unsigned int i) {EndType = i;}
	inline void setID(unsigned int i) {SegmentID = i;}
	inline void setVesselID(unsigned int i) {VesselID = i;}
	inline unsigned int getVesselID() {return(VesselID);}
	
		
private:
    unsigned int SegmentID;
    unsigned int VesselID;
    Vect2 mu;
    double q;
    double a1;
    double a2;
    Vect2 R1;
    Vect2 R2;
    double e1;
    float ForegdIntensity;
    float BackgdIntensity;
    float MADEstimate;
    float Likelihood;
    unsigned int EndType;
    Vect2 last_pos;
	Vect2 last_dir;
	int SelfIntersectCounter;
    //std::string Comment;
	double MinimumWidth;
	double MaximumWidth;
	double last_min_width;
	double PROP;
	double PredictedStep;
	double AspectRatio;
	double StepSize;
};

#endif
/*EndType
0 = Not an end
11 = End of Image 	(reduced)
12 = End of Image	(absolute)
2 = Hitting other vessel
31 = Hitting itself
32 = Selfintersecting
33 = flipped more than 3 times
4 =	
5 =	
6 =
7 = Width Discontinuity
8 = Exceeded MIN and MAX vessel width
*/
