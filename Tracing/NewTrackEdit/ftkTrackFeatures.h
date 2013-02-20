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
		fprintf(fp, " %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f \n",scalars[DISTANCE_TO_1],scalars[DISTANCE],scalars[INST_SPEED],scalars[ANGLE_REL_TO_1],scalars[CHANGE_DISTANCE_TO_1],scalars[HAS_CONTACT_TO_2],scalars[DISPLACEMENT_VEC_X],scalars[DISPLACEMENT_VEC_Y],scalars[DISPLACEMENT_VEC_Z]);
	}
};
class TrackFeatures{
	public:
		std::vector<ftk::LabelImageFeatures> intrinsic_features;
		std::vector<ftk::TrackPointFeatures> tfeatures;
		enum{ AVG_SPEED, AVG_DIST_TO_1, AVG_ANGLE_REL_TO_1,CHANGE_DISTANCE_TO_1,CONTACT_TO_2, DISPLACEMENT_VEC_X, DISPLACEMENT_VEC_Y, DISPLACEMENT_VEC_Z, PATHLENGTH, CONFINEMENT_RATIO};
		float scalars[CONFINEMENT_RATIO+1];
		void Fprintf(FILE* fp1 = stdout,FILE *fp2 = stdout)
		{
		//	fprintf(fp,"TrackFeatures at each time point:\n");
			if(intrinsic_features.size()>2)
			{
			for(int counter=0; counter< ((intrinsic_features.size()<tfeatures.size())?(intrinsic_features.size()):(tfeatures.size())); counter++)
			{
				intrinsic_features[counter].Fprintf(fp1,3);
				tfeatures[counter].Fprintf(fp1);
			}

//			for(int counter=0; counter<tfeatures.size(); counter++)
//			{
//				tfeatures[counter].Fprintf(fp1);
//			}
			fprintf(fp2,"%d %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n",intrinsic_features[0].num,\
					scalars[AVG_SPEED],scalars[AVG_DIST_TO_1],scalars[AVG_ANGLE_REL_TO_1],scalars[CHANGE_DISTANCE_TO_1],scalars[CONTACT_TO_2],\
					scalars[DISPLACEMENT_VEC_X],scalars[DISPLACEMENT_VEC_Y],scalars[DISPLACEMENT_VEC_Z],scalars[PATHLENGTH],\
					1/scalars[CONFINEMENT_RATIO]);
			}
		}
};

};
#endif
