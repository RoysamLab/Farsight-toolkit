/** @file SeedContainer3D.h
*   @brief Basic class for detecting seed points in 3D volume.
*
*   @author Alex James Tyrrell, Amit Mukherjee,
*/
#ifndef SEED_CONTAINER3D_H
#define SEED_CONTAINER3D_H

#include "itkImage.h"

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "TraceConfig.h"
#include "SeedPoint3D.h"
#include <algorithm>

typedef itk::Image<float, 3> ImageType3D;
typedef std::vector <SeedPoint3D *> SeedContainter3DType;

class SeedContainer3D: public itk::LightObject {

public:
	typedef SeedContainer3D Self;
	typedef  itk::SmartPointer< Self >  Pointer;
	typedef  itk::SmartPointer< const Self >  ConstPointer;
	itkNewMacro(Self);


	SeedPoint3D* getSeed(unsigned int i) {return(SeedContainer[i]);}
	unsigned int getNumberOfSeeds() {return(SeedContainer.size());}

	void Configure(TraceConfig::Pointer);
	void Detect(ImageType3D::Pointer, ImageType2D::Pointer );
	SeedContainter3DType getContainer() {return SeedContainer;}

private:
	long GridSpacing;
	int IntensityThreshold;
	SeedContainter3DType SeedContainer;

	void AddSeedPoint(ImageType2D::RegionType::IndexType, const float&);
	void AddSeedPoint(ImageType3D::IndexType , ImageType3D::PixelType , long );
	void Detect2Dseeds(ImageType2D::Pointer);
	void Detect3Dseeds(ImageType3D::Pointer);
	void MIPGenerator(ImageType3D::Pointer, ImageType2D::Pointer&);
	void LocateZPosition(ImageType3D::Pointer );



protected:
	SeedContainer3D(); //copy the config informations
	~SeedContainer3D();


};

#endif
