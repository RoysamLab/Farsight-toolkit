#include "ftkTrackFeatures.h"

namespace ftk
{
TrackPointFeatureInfoType TrackPointFeatures::Info[M]={
		{ "distance_to_1", "units", "description" },
		{ "distance", "units", "description" },
		{ "inst_speed", "units", "description" },
		{ "angle_rel_to_1", "units", "description" },
		{ "change_distance_to_1", "units", "description" },
		{ "has_contact_to_2", "units", "description" },
		{ "displacement_vec_x", "units", "description" },
		{ "displacement_vec_y", "units", "description" },
		{ "displacement_vec_z", "units", "description" },
};
TimeFeatureInfoType TrackFeatures::TimeInfo[NF]={
		{ "avg_speed", "units", "description" },
		{ "max_speed", "units", "description" },
		{ "min_speed", "units", "description" },
		{ "avg_dis_to_1", "units", "description" },
		{ "avg_angle_rel_to_1", "units", "description" },
		{ "change_distance_to_1", "units", "description" },
		{ "contact_to_2", "units", "description" },
		{ "displacement_vec_x", "units", "description" },
		{ "displacement_vec_y", "units", "description" },
		{ "displacement_vec_z", "units", "description" },
		{ "pathlength", "units", "description" },
		{ "total_distance", "units", "description" },
		{ "confinement_ratio", "units", "description" },
};

void TrackPointFeatures::Fprintf(FILE *fp)
{
	fprintf(fp, " %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f \n",scalars[DISTANCE_TO_1],scalars[DISTANCE],scalars[INST_SPEED],scalars[ANGLE_REL_TO_1],scalars[CHANGE_DISTANCE_TO_1],scalars[HAS_CONTACT_TO_2],scalars[DISPLACEMENT_VEC_X],scalars[DISPLACEMENT_VEC_Y],scalars[DISPLACEMENT_VEC_Z]);
}
void TrackFeatures::Fprintf(FILE* fp1,FILE *fp2)
{
	//if(intrinsic_features.size()>0)
	//{
	//	for(unsigned int counter=0; counter< ((intrinsic_features.size()<tfeatures.size())?(intrinsic_features.size()):(tfeatures.size())); counter++)
	//	{
	//		fprintf(fp1,"%d %d %0.3f %0.3f %0.3f ",intrinsic_features[counter].num, intrinsic_features[counter].time,intrinsic_features[counter].Centroid[0],intrinsic_features[counter].Centroid[1],intrinsic_features[counter].Centroid[2]);
	//		for(int counter1=0; counter1<= ftk::IntrinsicFeatures::CLUSTER_PROMINENCE; counter1++)
	//		{
	//			fprintf(fp1,"%0.3f ", intrinsic_features[counter].ScalarFeatures[counter1]);
	//		}
	//		tfeatures[counter].Fprintf(fp1);
	//	}

	//	fprintf(fp2,"%d %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n",intrinsic_features[0].num,\
	//			scalars[AVG_SPEED],scalars[AVG_DIST_TO_1],scalars[AVG_ANGLE_REL_TO_1],scalars[CHANGE_DISTANCE_TO_1],scalars[CONTACT_TO_2],\
	//			scalars[DISPLACEMENT_VEC_X],scalars[DISPLACEMENT_VEC_Y],scalars[DISPLACEMENT_VEC_Z],scalars[PATHLENGTH],\
	//			1/scalars[CONFINEMENT_RATIO]);
	//}
}
}