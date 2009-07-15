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

#ifndef __FTK_TRACK_FEATURES_H
#define __FTK_TRACK_FEATURES_H
#include <vector>
#include <stdio.h>
#include "ftkLabelImageToFeatures.h"

namespace ftk
{
class TrackPointFeatures{
	public:
	enum{DISTANCE_TO_1, DISTANCE, INST_SPEED, ANGLE_REL_TO_1, CHANGE_DISTANCE_TO_1, HAS_CONTACT_TO_2, DISPLACEMENT_VEC_X, DISPLACEMENT_VEC_Y, DISPLACEMENT_VEC_Z};
	float scalars[DISPLACEMENT_VEC_Z+1];
	void Fprintf(FILE *fp = stdout)
	{
		fprintf(fp, "%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n",scalars[DISTANCE_TO_1],scalars[DISTANCE],scalars[INST_SPEED],scalars[ANGLE_REL_TO_1],scalars[CHANGE_DISTANCE_TO_1],scalars[HAS_CONTACT_TO_2],scalars[DISPLACEMENT_VEC_X],scalars[DISPLACEMENT_VEC_Y],scalars[DISPLACEMENT_VEC_Z]);
	}
};
class TrackFeatures{
	public:
		std::vector<ftk::LabelImageFeatures> intrinsic_features;
		std::vector<ftk::TrackPointFeatures> tfeatures;
		enum{ AVG_SPEED, AVG_DIST_TO_1, AVG_ANGLE_REL_TO_1,CHANGE_DISTANCE_TO_1,CONTACT_TO_2, DISPLACEMENT_VEC_X, DISPLACEMENT_VEC_Y, DISPLACEMENT_VEC_Z, PATHLENGTH, CONFINEMENT_RATIO};
		float scalars[CONFINEMENT_RATIO+1];
		void Fprintf(FILE* fp = stdout)
		{
			fprintf(fp,"TrackFeatures at each time point:\n");
			for(int counter=0; counter< intrinsic_features.size(); counter++)
			{
				intrinsic_features[counter].Fprintf(fp);
			}
			for(int counter=0; counter<tfeatures.size(); counter++)
			{
				tfeatures[counter].Fprintf(fp);
			}
			fprintf(fp,"Overall Track Features : %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n",\
					scalars[AVG_SPEED],scalars[AVG_DIST_TO_1],scalars[AVG_ANGLE_REL_TO_1],scalars[CHANGE_DISTANCE_TO_1],scalars[CONTACT_TO_2],\
					scalars[DISPLACEMENT_VEC_X],scalars[DISPLACEMENT_VEC_Y],scalars[DISPLACEMENT_VEC_Z],scalars[PATHLENGTH],\
					scalars[CONFINEMENT_RATIO]);
			printf("\n\n");
		}
};

};
#endif
