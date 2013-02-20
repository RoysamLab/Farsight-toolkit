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

SubTrackFeatureInfoType SubTrackFeatures::SubTrackInfo[NF]={
		{ "dis_x", "units", "description" },
		{ "dis_y", "units", "description" },
		{ "dis_z", "units", "description" },
		{ "length", "units", "description" },
		{ "cos_angle_rel_prev_dir", "units", "description" },
		{ "cos_angle_rel_path_dir", "units", "description" },
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
void TrackingFeatureComputation::ComputeSubTrackFeatures(std::vector< std::vector <ftk::TrackPoint> > * track_points)
{
	printf("Started ComputeSubTrackFeatures ...\n");
	typedef ftk::SubTrackFeatures STF;
	unsigned int n_tracks = track_points->size();
	std::vector<double> disx_vec;
	std::vector<double> disy_vec;
	std::vector<double> disz_vec;
	std::vector<double> length_vec;

	FILE * fp1 = fopen("C:\\Lab\\AminFiles\\Debug\\SubTrackFeaturesTests1.txt","w");

	for(unsigned int i=0; i<n_tracks; ++i)			// iterate over number of tracks
	{
		std::vector<ftk::TrackPoint> * current_track = &(*track_points)[i];
	
		if(current_track->size()<=1)
			continue;

		std::vector< STF > shared_id_tracks;
		unsigned int n_time_points = current_track->size();
		for(unsigned int j=0; j<n_time_points-1; ++j)		// iterate over time
		{
			ftk::TrackPoint * curr_bit = &(*current_track)[j];
			ftk::TrackPoint * next_bit = &(*current_track)[j+1];
			STF stf;
			stf.num = curr_bit->id;
			stf.t1 = curr_bit->t;
			stf.t2 = next_bit->t;

			float timediff = (float)(next_bit->t - curr_bit->t);
			if(timediff > 0)																	// Later scale this thing with the spacing
			{
				timediff = 1/timediff;

				// distance computations:
				//stf.scalars[STF::DIS_X] = (float)(next_bit->x - curr_bit->x)*timediff;
				//stf.scalars[STF::DIS_Y] = (float)(next_bit->y - curr_bit->y)*timediff;
				//stf.scalars[STF::DIS_Z] = (float)(next_bit->z - curr_bit->z)*timediff;
				stf.scalars[STF::DIS_X] = (float)(next_bit->x - curr_bit->x);
				stf.scalars[STF::DIS_Y] = (float)(next_bit->y - curr_bit->y);
				stf.scalars[STF::DIS_Z] = (float)(next_bit->z - curr_bit->z);
				



				fprintf(fp1,"%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",curr_bit->id,stf.t1,stf.t2,stf.scalars[STF::DIS_X],curr_bit->BoundingBox[0], curr_bit->BoundingBox[1],\
																					curr_bit->BoundingBox[2], curr_bit->BoundingBox[3],stf.scalars[STF::DIS_Y],next_bit->BoundingBox[0],next_bit->BoundingBox[1],\
																					next_bit->BoundingBox[2],next_bit->BoundingBox[3]);

		
				disx_vec.push_back((double)stf.scalars[STF::DIS_X]);									// store the data to compute the mean and standard deviation
				disy_vec.push_back((double)stf.scalars[STF::DIS_Y]);
				disz_vec.push_back((double)stf.scalars[STF::DIS_Z]);

				// length computations:
				stf.scalars[STF::LENGTH] = sqrt(stf.scalars[STF::DIS_X]*stf.scalars[STF::DIS_X] + stf.scalars[STF::DIS_Y]*stf.scalars[STF::DIS_Y] + stf.scalars[STF::DIS_Z]*stf.scalars[STF::DIS_Z]);
				length_vec.push_back(stf.scalars[STF::LENGTH]);
				
				shared_id_tracks.push_back(stf);
			}
			else
			{
				printf("Fatal Error: timediff is %d",timediff);
				scanf("%*d");
			}
		}

		_sub_track_feats.push_back(shared_id_tracks);
	}
	fclose(fp1);
	for(unsigned int i=0; i<_sub_track_feats.size(); ++i)
	{
		float tot_disp_X = 0.0;
		float tot_disp_Y = 0.0;
		float tot_disp_Z = 0.0;

		tot_disp_X += _sub_track_feats[i][0].scalars[STF::DIS_X];
		tot_disp_Y += _sub_track_feats[i][0].scalars[STF::DIS_X];
		tot_disp_Z += _sub_track_feats[i][0].scalars[STF::DIS_X];

		for(unsigned int j=1; j<_sub_track_feats[i].size(); ++j)
		{
			_sub_track_feats[i][j].scalars[STF::COS_ANGLE_REL_TO_PREV] = _sub_track_feats[i][j-1].scalars[STF::DIS_X]*_sub_track_feats[i][j].scalars[STF::DIS_X]+\
																		 _sub_track_feats[i][j-1].scalars[STF::DIS_Y]*_sub_track_feats[i][j].scalars[STF::DIS_Y]+\
																		 _sub_track_feats[i][j-1].scalars[STF::DIS_Z]*_sub_track_feats[i][j].scalars[STF::DIS_Z];

			_sub_track_feats[i][j].scalars[STF::COS_ANGLE_REL_TO_PREV]/=_sub_track_feats[i][j-1].scalars[STF::LENGTH];
			_sub_track_feats[i][j].scalars[STF::COS_ANGLE_REL_TO_PREV]/=_sub_track_feats[i][j].scalars[STF::LENGTH];
			tot_disp_X +=  _sub_track_feats[i][j].scalars[STF::DIS_X];
			tot_disp_Y +=  _sub_track_feats[i][j].scalars[STF::DIS_Y];
			tot_disp_Z +=  _sub_track_feats[i][j].scalars[STF::DIS_Z];

		}
		float tot_disp_Length = sqrt((tot_disp_X*tot_disp_X)+(tot_disp_Y*tot_disp_Y)+(tot_disp_Z*tot_disp_Z));
		if(tot_disp_Length==0.0)
		{
			tot_disp_Length = 1.0;		// FIX ME
		}
		for(unsigned int j=0; j<_sub_track_feats[i].size(); ++j)
		{
			_sub_track_feats[i][j].scalars[STF::COS_ANGLE_REL_PATH_DIR] = _sub_track_feats[i][j].scalars[STF::DIS_X]*tot_disp_X+\
																	      _sub_track_feats[i][j].scalars[STF::DIS_Y]*tot_disp_Y+\
																	      _sub_track_feats[i][j].scalars[STF::DIS_Z]*tot_disp_Z;
			_sub_track_feats[i][j].scalars[STF::COS_ANGLE_REL_PATH_DIR] /=_sub_track_feats[i][j].scalars[STF::LENGTH];
			_sub_track_feats[i][j].scalars[STF::COS_ANGLE_REL_PATH_DIR] /=tot_disp_Length;
		}



	}

	//FILE * fp = fopen("C:\\Lab\\AminFiles\\Debug\\SubTrackFeaturesTests.txt","w");
	//for(unsigned int i=0; i<_sub_track_feats.size(); ++i)
	//{
	//	for(unsigned int j=0; j<_sub_track_feats[i].size(); ++j)
	//	{
	//		fprintf(fp,"%d\t%d\t",_sub_track_feats[i][j].num,_sub_track_feats[i][j].t1);
	//		printf("%d\t%d\t",_sub_track_feats[i][j].num,_sub_track_feats[i][j].t1);
	//		for(int k=0; k<STF::COS_ANGLE_REL_PATH_DIR+1; ++k)
	//		{
	//			fprintf(fp,"%f\t",_sub_track_feats[i][j].scalars[k]);
	//		}
	//		fprintf(fp,"\n");
	//	}
	//}
	//fclose(fp);
	printf("Finished ComputeSubTrackFeatures ...\n");



}

}