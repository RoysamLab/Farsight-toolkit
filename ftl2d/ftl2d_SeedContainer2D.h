#ifndef SEED_CONTAINER2D_H
#define SEED_CONTAINER2D_H


#include "ftl2d_TraceConfig.h"
#include "ftl2d_SeedPoint2D.h"
#include <algorithm>


class SeedContainer2D {
	public:
	SeedContainer2D( TraceConfig * ); //copy the config informations
	~SeedContainer2D();
	void AddSeedPoint(ImageType::RegionType::IndexType ndx, const float&);
	SeedPoint2D* getSeed(unsigned int i) {return(SeedContainer[i]);}
	unsigned int getNumberOfSeeds() {return(SeedContainer.size());}

	void Detect(ImageType::Pointer im);
	void Validate(ImageType::Pointer im);
	void SortSeeds();
//	bool SortIntensityFunction(const SeedPoint2D* &, const SeedPoint2D* & );

	private:
	int GridSpacing;
	int IntensityThreshold;
	unsigned char valid;
	std::vector <SeedPoint2D *> SeedContainer;
};

#endif
