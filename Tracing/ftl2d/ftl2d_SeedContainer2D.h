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
