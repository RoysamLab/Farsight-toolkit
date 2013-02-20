#ifndef __FTK_TRACK_FEATURES_H
#define __FTK_TRACK_FEATURES_H
#include <vector>
#include <stdio.h>
#include <math.h>

#include "ftkLabelImageToFeatures.h"
#include <ftkFeatures/ftkObject.h>
//#include "ftkIntrinsicFeatures.h"
#define MAX(a,b) (((a) > (b))?(a):(b))


namespace ftk
{
typedef struct {
	std::string name;
	std::string units;
	std::string description;
} TrackPointFeatureInfoType;
typedef struct {
	std::string name;
	std::string units;
	std::string description;
} TimeFeatureInfoType;
typedef struct {
	std::string name;
	std::string units;
	std::string description;
} SubTrackFeatureInfoType;

typedef struct{
	double x,y,z;
	int id,t;
	float BoundingBox[6];
} TrackPoint;





class TrackPointFeatures{
public:
	TrackPointFeatures(){
		for(int counter=0; counter <DISPLACEMENT_VEC_Z+1; counter++)
		{
			scalars[counter] = 0.0;
		}
	}
	enum{DISTANCE_TO_1, DISTANCE, INST_SPEED, ANGLE_REL_TO_1, CHANGE_DISTANCE_TO_1, HAS_CONTACT_TO_2, DISPLACEMENT_VEC_X, DISPLACEMENT_VEC_Y, DISPLACEMENT_VEC_Z};

	static const int M = DISPLACEMENT_VEC_Z+1;	//This is the number of scalar track features 
	float scalars[M];
	int num;

	static TrackPointFeatureInfoType Info[M];
	void Fprintf(FILE *fp = stdout);

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

	
		std::vector<ftk::IntrinsicFeatures> intrinsic_features;// figure out what makes some of the vectors empty.
		std::vector<ftk::TrackPointFeatures> tfeatures;
		enum{ AVG_SPEED, MAX_SPEED, MIN_SPEED, AVG_DIST_TO_1, AVG_ANGLE_REL_TO_1,CHANGE_DISTANCE_TO_1, CONTACT_TO_2, DISPLACEMENT_VEC_X,\
			  DISPLACEMENT_VEC_Y, DISPLACEMENT_VEC_Z, PATHLENGTH, TOTAL_DISTANCE, CONFINEMENT_RATIO};
		static const int NF = CONFINEMENT_RATIO+1;
		float scalars[NF];
		static TimeFeatureInfoType TimeInfo[NF];

		void Fprintf(FILE* fp1 = stdout,FILE *fp2 = stdout);

};

class SubTrackFeatures
{
	public:
		SubTrackFeatures()
		{
			for(int counter=0; counter< COS_ANGLE_REL_PATH_DIR+1; counter++)
			{
				scalars[counter] = 0.0;
			}
		}
		enum{ DIS_X, DIS_Y, DIS_Z, LENGTH, COS_ANGLE_REL_TO_PREV, COS_ANGLE_REL_PATH_DIR};

		static const int NF = COS_ANGLE_REL_PATH_DIR+1;
		float scalars[NF];
		static SubTrackFeatureInfoType SubTrackInfo[NF];
		int num;
		int t1,t2;

			
};

class TrackingFeatureComputation
{	
	public:
		TrackingFeatureComputation()
		{
		}

		
		void ComputeSubTrackFeatures(std::vector< std::vector <ftk::TrackPoint> > * track_points);
		std::vector< std::vector<ftk::SubTrackFeatures> > GetSubTrackFeatures(void){return this->_sub_track_feats;};

	private:
		std::vector< std::vector< ftk::SubTrackFeatures > > _sub_track_feats;


};

}
#endif
