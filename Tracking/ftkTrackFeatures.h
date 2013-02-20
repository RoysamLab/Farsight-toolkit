#ifndef __FTK_TRACK_FEATURES_H
#define __FTK_TRACK_FEATURES_H
#include <vector>
#include <stdio.h>
#include "ftkLabelImageToFeatures.h"
#include "ftkIntrinsicFeatures.h"



namespace ftk
{
class TrackPointFeatures{
public:
	TrackPointFeatures(){
		for(int counter=0; counter <DISPLACEMENT_VEC_Z+1; counter++)
		{
			scalars[counter] = 0.0;
		}
	}
	enum{DISTANCE_TO_1, DISTANCE, INST_SPEED, ANGLE_REL_TO_1, CHANGE_DISTANCE_TO_1, HAS_CONTACT_TO_2, DISPLACEMENT_VEC_X, DISPLACEMENT_VEC_Y, DISPLACEMENT_VEC_Z};
	float scalars[DISPLACEMENT_VEC_Z+1];
	void Fprintf(FILE *fp = stdout)
	{
		fprintf(fp, " %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f \n",scalars[DISTANCE_TO_1],scalars[DISTANCE],scalars[INST_SPEED],scalars[ANGLE_REL_TO_1],scalars[CHANGE_DISTANCE_TO_1],scalars[HAS_CONTACT_TO_2],scalars[DISPLACEMENT_VEC_X],scalars[DISPLACEMENT_VEC_Y],scalars[DISPLACEMENT_VEC_Z]);
	}
};
class TrackFeatures{
	public:
		TrackFeatures()
		{
			for(int counter=0; counter< CONFINEMENT_RATIO+1; counter++)
			{
				scalars[counter] = 0.0;
			}
		}
	
		std::vector<ftk::IntrinsicFeatures> intrinsic_features;
		std::vector<ftk::TrackPointFeatures> tfeatures;
		enum{ AVG_SPEED, AVG_DIST_TO_1, AVG_ANGLE_REL_TO_1,CHANGE_DISTANCE_TO_1, CONTACT_TO_2, DISPLACEMENT_VEC_X, DISPLACEMENT_VEC_Y, DISPLACEMENT_VEC_Z, PATHLENGTH, CONFINEMENT_RATIO};
		float scalars[CONFINEMENT_RATIO+1];
		void Fprintf(FILE* fp1 = stdout,FILE *fp2 = stdout)
		{
		//	fprintf(fp,"TrackFeatures at each time point:\n");
			if(intrinsic_features.size()>0)
			{
			for(unsigned int counter=0; counter< ((intrinsic_features.size()<tfeatures.size())?(intrinsic_features.size()):(tfeatures.size())); counter++)
			{
				fprintf(fp1,"%d %d %0.3f %0.3f %0.3f ",intrinsic_features[counter].num, intrinsic_features[counter].time,intrinsic_features[counter].Centroid[0],intrinsic_features[counter].Centroid[1],intrinsic_features[counter].Centroid[2]);
				for(int counter1=0; counter1<= ftk::IntrinsicFeatures::CLUSTER_PROMINENCE; counter1++)
				{
					fprintf(fp1,"%0.3f ", intrinsic_features[counter].ScalarFeatures[counter1]);
				}
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
}
#endif
