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

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef SNAKEOPERATION_H
#define SNAKEOPERATION_H

#include "ImageOperation.h"
#include "PointOperation.h"
#include "itkLinearInterpolateImageFunction.h"

bool isinf(float x);
bool isnan(float x);
float norm_density(float x, float mu, float sigma);

class SnakeListClass;

class SnakeClass 
{
public:

    SnakeClass(void);
	SnakeClass operator=(SnakeClass Snake);

	PointList3D Cu;
	PointList3D Cu_Backup;
	PointList3D BranchPt;
	PointList3D RootPt;
	
	//vnl_vector<float> Ru;
	std::vector<float> Ru;

	//PointList3D Probe;

	Point3D temp_pt;
	Point3D head_pt;
	Point3D tail_pt;

	int collision;
	int tail_collision_snake_id;
	int head_collision_snake_id;
    bool hit_boundary;

	ImageOperation *IM;
    
	SnakeListClass *SnakeList;

	void Branch_Adjustment();
	void Nail_Branch();

	bool Jump_Over_Gap(int dist, int angle, int step, int ratio, int ht);
	bool Jump_Over_Crossover(int dist, int angle, int step, int ratio, int collision);
	bool Jump_Over_Crossover_New(int dist, int collision);

	void SetTracedSnakes(SnakeListClass *S);
	void SetImage(ImageOperation *I_Input);
	void Set_Seed_Point(PointList3D seeds);
    void Set_Seed_Point(Point3D seed);

	void Skeleton_Expansion();
	void Expand_Seed_Point(int expand_distance);

    void Grow_Snake_Point();
	void OpenSnakeDeform(float alpha, int ITER, int pt_distance, float beta, float kappa, float gamma, int N_Active, bool freeze_body);
	void OpenSnakeStretch(float alpha, int ITER, int pt_distance, float beta, float kappa, float gamma, 
		    float stretchingRatio, int collision_dist, int N_Active, int minimum_length, bool automatic_merging, 
			int max_angle, bool freeze_body, int s_force, int snake_id);
   
	void OpenSnakeStretch_4D(float alpha, int ITER, int pt_distance, float beta, float kappa, float gamma, 
                                float stretchingRatio, int collision_dist, int minimum_length, 
								bool automatic_merging, int max_angle, bool freeze_body, int s_force, int snake_id, int tracing_model);

	bool Check_Validity(float minimum_length, float repeat_ratio, int repeat_dist, int snake_id);
	bool Check_Head_Collision(ImageType::IndexType in, int collision_dist, int minimum_length, bool automatic_merging, int max_angle, int snake_id);
	bool Check_Tail_Collision(ImageType::IndexType in, int collision_dist, int minimum_length, bool automatic_merging, int max_angle, int snake_id);
	bool Check_Head_Leakage(ImageType::IndexType in);
    bool Check_Tail_Leakage(ImageType::IndexType in);
	bool Compute_Seed_Force(int head_tail, int distance);
	void Estimate_Radius();

	vnl_matrix<float> makeOpenA(float alpha, float beta, int N);
};

class SnakeListClass
{
public:

    SnakeListClass(void);
	SnakeListClass operator= (SnakeListClass SnakeList);
    
    int NSnakes;
	SnakeClass *Snakes;
	vnl_vector<int> valid_list;
    
	ImageOperation *IM;

	void AddSnake(SnakeClass snake);
	void RemoveSnake(int idx);
	void SplitSnake(int idx_snake, int idx_pt);
	void MergeSnake(int idx1, int idx2, bool im_coding);
	void CreateBranch(int idx1, int idx2);
	SnakeClass GetSnake(int idx);
	void SetNSpace(int N);
	void SetImage(ImageOperation *I_Input);

};

struct SnakeTree
{
	Point3D root_point;
    PointList3D branch_point;
	std::vector<int> parent_list;
	PointList3D points;
	std::vector<float> Ru;
	std::vector<int> snake_id;
};

#endif
