#ifndef SEED2SEG_H_
#define SEED2SEG_H_

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "itkImage.h"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <vector>
#include <algorithm>
#include <cstdlib>
#include "vnl/vnl_vector_fixed.h"


#include "TVessel.h"
#include "SegInit.h"
#include "SeedContainer3D.h"

typedef itk::Image<float , 3> ImageType3D;
typedef std::vector<SeedPoint3D *> SeedContainerType;
typedef std::vector <TVessel *> StartSegContainerType;

class Seed2Seg: public itk::LightObject	{
public:
	typedef Seed2Seg Self;
	typedef  itk::SmartPointer< Self >  Pointer;
	typedef  itk::SmartPointer< const Self >  ConstPointer;
	itkNewMacro(Self);
	typedef vnl_vector_fixed<double,3> Vect3;

	void ComuputeStartSegments(SeedContainer3D::Pointer , ImageType3D::Pointer, TraceConfig::Pointer);
	void SortStartSegments();
	unsigned int getNumberOfStartSegments() {return SSContainer.size();}
	TVessel* getStartSegment(unsigned int i) {return SSContainer[i];}
	
private:
	bool HitTestSeedPoints(SeedPoint3D*);
	bool HitTestSeg(TVessel*  );
	//bool Sortfunction(TVessel*s1, TVessel*s2);

	double iterations;
	double AS_RATIO;
	StartSegContainerType SSContainer;



protected:
	Seed2Seg();

	~Seed2Seg();
};

#endif /*SEGINIT_H_*/
