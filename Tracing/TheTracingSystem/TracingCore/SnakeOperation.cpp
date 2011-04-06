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
#include "SnakeOperation.h"

SnakeClass::SnakeClass(void)
{
	//Cu.SetN(1024);
	collision = 0;
	hit_boundary = false;
	tail_collision_snake_id = -1;
	head_collision_snake_id = -1;
}

SnakeClass SnakeClass::operator=(SnakeClass Snake)
{ 
   //Cu.SetN(1024);
   Cu = Snake.Cu;
   BranchPt = Snake.BranchPt;
   RootPt = Snake.RootPt;
   Ru = Snake.Ru;
   
   return *this;
}

void SnakeClass::SetImage(ImageOperation *I_Input)
{
   IM = I_Input;
}

void SnakeClass::Set_Seed_Point(PointList3D seeds)
{
   collision = 0;
   hit_boundary = false;
   //Cu.NP = 0;
   Cu.RemoveAllPts();
   Cu.AddPtList(seeds);
}

void SnakeClass::Set_Seed_Point(Point3D seed)
{
   collision = 0;
   hit_boundary = false;
   //Cu.NP = 0;
   Cu.RemoveAllPts();
   Cu.AddPt(seed);
}

void SnakeClass::Nail_Branch()
{
  //replace head or tail with new branch point
   for( int i = 0; i < BranchPt.GetSize(); i++ )
   {
	   BranchPt.Pt[i].check_out_of_range_3D(IM->SM,IM->SN,IM->SZ);
	   if( BranchPt.Pt[i].GetDistTo(Cu.GetFirstPt()) > BranchPt.Pt[i].GetDistTo(Cu.GetLastPt()) ) 
	   {
	      Cu.Pt[Cu.NP-1] = BranchPt.Pt[i];
	   }
	   else
	   {
	      Cu.Pt[0] = BranchPt.Pt[i];
	   }
   }
}

void SnakeClass::Branch_Adjustment()
{

   int iter_num = 10;
   typedef itk::VectorLinearInterpolateImageFunction< 
                       GradientImageType, float >  GradientInterpolatorType;

   GradientInterpolatorType::Pointer interpolator = GradientInterpolatorType::New();
   interpolator->SetInputImage(IM->IGVF);

   for( int i = 0; i < BranchPt.GetSize(); i++ )
   {
	 int j = 0;
     while( j < iter_num )
	 {
	    GradientImageType::IndexType index; 
	    BranchPt.Pt[i].check_out_of_range_3D(IM->SM,IM->SN,IM->SZ);
		index[0] = (BranchPt.Pt[i].x);
		index[1] = (BranchPt.Pt[i].y);
		index[2] = (BranchPt.Pt[i].z);
	    GradientPixelType gradient = interpolator->EvaluateAtIndex(index);
		temp_pt.x = gradient[0];
        temp_pt.y = gradient[1];
		temp_pt.z = gradient[2];
        BranchPt.Pt[i] = BranchPt.Pt[i] + temp_pt;
		j++;
	 }
   }

  Nail_Branch();
 
}

void SnakeClass::Skeleton_Expansion()
{
   if( Cu.NP != 1 )
	   return;
	   
   int SM = IM->I->GetLargestPossibleRegion().GetSize()[0];
   int SN = IM->I->GetLargestPossibleRegion().GetSize()[1];
   int SZ = IM->I->GetLargestPossibleRegion().GetSize()[2];

   typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
   typedef itk::NeighborhoodIterator< ImageType > NeighborhoodIteratorType;
   typedef itk::NeighborhoodIterator< LabelImageType > NeighborhoodIteratorType_Label;


   NeighborhoodIteratorType::RadiusType radius;
   NeighborhoodIteratorType_Label::RadiusType radius_label;
   radius.Fill(1);
   radius_label.Fill(1);
   NeighborhoodIteratorType it(radius, IM->VBW, IM->VBW->GetRequestedRegion());
   NeighborhoodIteratorType_Label it_label(radius_label, IM->IL, IM->IL->GetRequestedRegion());
   NeighborhoodIteratorType_Label it_label_t(radius_label, IM->IL_Tracing, IM->IL_Tracing->GetRequestedRegion());
 
   ImageType::IndexType index;
   index[0] = Cu.Pt[0].x;
   index[1] = Cu.Pt[0].y;
   index[2] = Cu.Pt[0].z;
   it.SetLocation(index);
   it_label.SetLocation(index);
   it_label_t.SetLocation(index);

   //if current seed point is not in skeleton, move it to nearby skeleton point
   NeighborhoodIteratorType::OffsetType nextMove;
   nextMove.Fill(0);

   if( it.GetCenterPixel() == 0 )
   {
      for( unsigned int i = 0; i < it.Size(); i++ )
	  {
	    if( it.GetPixel(i) == 1 && it_label.GetPixel(i) == 0 )
		{
		  index[0] = it.GetIndex(i)[0];
		  index[1] = it.GetIndex(i)[1];
		  index[2] = it.GetIndex(i)[2];
		  it.SetLocation(index);
		  it_label.SetLocation(index);
          it_label_t.SetLocation(index);
		  break;
 		}
	  }
   }

   ImageType::IndexType start_idx = it.GetIndex(13);
   //check the validity of start_idx
   if( start_idx[0] < 0 || start_idx[1] < 0 || start_idx[2] < 0 || 
	   start_idx[0] >= SM || start_idx[1] >= SN || start_idx[2] >= SZ )
   {
     Expand_Seed_Point(3);
	 return;
   }

   //if cannot find nearby skeleton point, use normal seed exapnsion instead
   if( it.GetCenterPixel() == 0 )
   {
      Expand_Seed_Point(3);
	  return;
   }

   //if current point has more than three neighbor skeleton points, use normal seed exapnsion instead
   int skeleton_neighbor = 0;
   for( unsigned int i = 0; i < it.Size(); i++ )
   {
        if( i == 13 )
			continue;

	    if( it.GetPixel(i) == 1 && it_label.GetPixel(i) == 0 )
		{
		  skeleton_neighbor++;
 		}
	}

   if( skeleton_neighbor >= 3 || skeleton_neighbor == 0 )
   {
      Expand_Seed_Point(3);
	  return;
   }

   skeleton_neighbor = 1;

   //trace along one direction first
   while( skeleton_neighbor == 1 )
   {
      nextMove.Fill(0);
      for( unsigned int i = 0; i < it.Size(); i++ )
	  {
	    if( i == 13 )
			continue;
	    if( it.GetPixel(i) == 1 && it_label.GetPixel(i) == 0 && it_label_t.GetPixel(i) != 1 )
		{
			nextMove = it.GetOffset(i);
			it_label_t.SetCenterPixel( 1 ); 

			//add points
			ImageType::IndexType idx = it.GetIndex(13);
			temp_pt.x = idx[0];
			temp_pt.y = idx[1];
			temp_pt.z = idx[2];
			Cu.AddPt(temp_pt);
			//Cu.Pt[Cu.NP].x = idx[0];
            //Cu.Pt[Cu.NP].y = idx[1];
			//Cu.Pt[Cu.NP].z = idx[2];
			if( Cu.Pt[Cu.NP-1].check_out_of_range_3D(SM,SN,SZ) )
			{ 
			  Cu.RemovePt();
		      return;
			}
			//else
			//{
            //  Cu.NP++;
			//}

			it += nextMove;
			it_label += nextMove;
			it_label_t += nextMove;	

			break;
		}
	  }

	  if( nextMove[0] == 0 && nextMove[1] == 0 && nextMove[2] == 0 )
	  {
	      return;
	  }

		//find skeleton neighbors for current index
       skeleton_neighbor = 0;
	   for( unsigned int i = 0; i < it.Size(); i++ )
	   {
	    if( i == 13 )
			continue;
	    if( it.GetPixel(i) == 1 && it_label.GetPixel(i) == 0 && it_label_t.GetPixel(i) != 1 )
		{
			skeleton_neighbor++;
		}
	   }
	   //add end point before tracing ends
	   if( skeleton_neighbor == 0)
	   {   
		    //add points
		    ImageType::IndexType idx = it.GetIndex(13);
			temp_pt.x = idx[0];
			temp_pt.y = idx[1];
			temp_pt.z = idx[2];
			Cu.AddPt(temp_pt);
			//Cu.Pt[Cu.NP].x = idx[0];
            //Cu.Pt[Cu.NP].y = idx[1];
			//Cu.Pt[Cu.NP].z = idx[2];
			if( Cu.Pt[Cu.NP-1].check_out_of_range_3D(SM,SN,SZ) )
			{ 
			  Cu.RemovePt();
		      return;
			}
			//else
			//{
            //  Cu.NP++;
			//}
	   }
   }

   //trace along the other (if exists) direction
   Cu.Flip();

   it.SetLocation(start_idx);
   it_label.SetLocation(start_idx);
   it_label_t.SetLocation(start_idx);

   skeleton_neighbor = 1;
   while( skeleton_neighbor == 1 )
   {
       nextMove.Fill(0);
      for( unsigned int i = 0; i < it.Size(); i++ )
	  {
	    if( i == 13 )
			continue;
	    if( it.GetPixel(i) == 1 && it_label.GetPixel(i) == 0 && it_label_t.GetPixel(i) != 1 )
		{
			nextMove = it.GetOffset(i);
			
			//add points
		   if( it_label_t.GetCenterPixel() != 1 )
		   {
			ImageType::IndexType idx = it.GetIndex(13);
			temp_pt.x = idx[0];
			temp_pt.y = idx[1];
			temp_pt.z = idx[2];
			Cu.AddPt(temp_pt);
			//Cu.Pt[Cu.NP].x = idx[0];
            //Cu.Pt[Cu.NP].y = idx[1];
			//Cu.Pt[Cu.NP].z = idx[2];
			if( Cu.Pt[Cu.NP-1].check_out_of_range_3D(SM,SN,SZ) )
			{ 
			  Cu.RemovePt();
		      return;
			}
			//else
			//{
            //  Cu.NP++;
			//}
		   }

			it_label_t.SetCenterPixel( 1 ); 

			it += nextMove;
			it_label += nextMove;
			it_label_t += nextMove;	

			break;
		}
	  }

	  if( nextMove[0] == 0 && nextMove[1] == 0 && nextMove[2] == 0 )
	  {
	      return;
	  }

	   //find skeleton neighbors for current index
       skeleton_neighbor = 0;
	   for( unsigned int i = 0; i < it.Size(); i++ )
	   {
	    if( i == 13 )
			continue;
	    if( it.GetPixel(i) == 1 && it_label.GetPixel(i) == 0 && it_label_t.GetPixel(i) != 1 )
		{
			skeleton_neighbor++;
		}
	   }
	   //add end point before tracing ends
	   if( skeleton_neighbor == 0)
	   {   
		    //add points
		    ImageType::IndexType idx = it.GetIndex(13);
			temp_pt.x = idx[0];
			temp_pt.y = idx[1];
			temp_pt.z = idx[2];
			Cu.AddPt(temp_pt);
			//Cu.Pt[Cu.NP].x = idx[0];
            //Cu.Pt[Cu.NP].y = idx[1];
			//Cu.Pt[Cu.NP].z = idx[2];
			if( Cu.Pt[Cu.NP-1].check_out_of_range_3D(SM,SN,SZ) )
			{ 
			  Cu.RemovePt();
		      return;
			}
			//else
			//{
            //  Cu.NP++;
			//}
	   }
   }

}

void SnakeClass::Expand_Seed_Point(int expand_distance)
{
   int SM = IM->V1->GetLargestPossibleRegion().GetSize()[0];
   int SN = IM->V1->GetLargestPossibleRegion().GetSize()[1];
   int SZ = IM->V1->GetLargestPossibleRegion().GetSize()[2];
   //int expand_distance = 5;
   GradientImageType::IndexType index; 

  if( Cu.NP == 1 )
  {
    Cu.Pt[Cu.NP-1].check_out_of_range_3D(SM,SN,SZ);
	index[0] = (Cu.Pt[Cu.NP-1].x);
    index[1] = (Cu.Pt[Cu.NP-1].y);
	index[2] = (Cu.Pt[Cu.NP-1].z);
  
    GradientPixelType first_p_direction = IM->V1->GetPixel(index);
	temp_pt.x = first_p_direction[0];
    temp_pt.y = first_p_direction[1];
	temp_pt.z = first_p_direction[2];
	temp_pt = Cu.Pt[Cu.NP-1] - temp_pt * 2;
	Cu.AddTailPt(temp_pt);

    temp_pt.x = first_p_direction[0];
    temp_pt.y = first_p_direction[1];
	temp_pt.z = first_p_direction[2];
	temp_pt = Cu.Pt[Cu.NP-1] + temp_pt * 2;
	Cu.AddPt(temp_pt);

     for( int i = 0; i < Cu.GetSize(); i++ )
	 {
		 Cu.Pt[i].check_out_of_range_3D(SM,SN,SZ);
	 }
  }

   /*if(Cu.NP == 1)
   {
     for(int i = 1; i < expand_distance; i++)
	 {
	   Cu.Pt[Cu.NP-1].check_out_of_range_3D(SM,SN,SZ);
	   index[0] = (Cu.Pt[Cu.NP-1].x);
       index[1] = (Cu.Pt[Cu.NP-1].y);
	   index[2] = (Cu.Pt[Cu.NP-1].z);
	   GradientPixelType first_p_direction = IM->V1->GetPixel(index);
	   Point3D temp(first_p_direction[0],first_p_direction[1],first_p_direction[2]);
	   temp = Cu.Pt[Cu.NP-1] + temp * 1;
       Cu.AddPt(temp);
	 }
    
	 for( int i = 0; i < Cu.GetSize(); i++ )
	 {
		 Cu.Pt[i].check_out_of_range_3D(SM,SN,SZ);
	 }
   }*/

  //for(int j = 0; j < Cu.NP; j++)
  //{
  //   std::cout<< Cu.Pt[j].x<<","<<Cu.Pt[j].y<<","<<Cu.Pt[j].z<<std::endl;
  //}
}

bool SnakeClass::Jump_Over_Gap(int dist, int angle, int step, int ratio, int ht)
{
	bool jump_gap = false;

	int SM = IM->SM;
    int SN = IM->SN;
    int SZ = IM->SZ;

	//Cu_Backup.NP = 0;
	Cu_Backup.RemoveAllPts();

	int pi = 3.14159265;
	float x,y,z,temp_x,temp_y,temp_z;
    GradientImageType::IndexType index1,index2;
	GradientPixelType v1,v2;

   if( ht == 0 )
   {
	index1[0] = Cu.Pt[0].x;
	index1[1] = Cu.Pt[0].y;
	index1[2] = Cu.Pt[0].z;
	v1 = IM->V1->GetPixel(index1);

	x = Cu.Pt[0].x - Cu.Pt[1].x;
    y = Cu.Pt[0].y - Cu.Pt[1].y;
	z = Cu.Pt[0].z - Cu.Pt[1].z;
    float theta = acos(z/sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

	float phi = atan2(y,x);

	int k = 0;

	int step_size = 2*angle/step;

	vnl_vector<float> weight(10000);
	weight.fill(0);

	for( int j = 1; j <= dist; j++ )
	{
      for( int i = -angle; i <= angle; i += step_size)
	  {
	    for( int z = -angle/ratio; z <= angle/ratio; z += (step_size/ratio))
		{
		 temp_x = j * sin(theta + (float)z * pi/180) * cos(phi + (float)i * pi/180);
	     temp_y = j * sin(theta + (float)z * pi/180) * sin(phi + (float)i * pi/180);
	     temp_z = j * cos(theta + (float)z * pi/180);
	     //Cu_Backup.Pt[k].x = Cu.Pt[0].x + temp_x;
         //Cu_Backup.Pt[k].y = Cu.Pt[0].y + temp_y;
	     //Cu_Backup.Pt[k].z = Cu.Pt[0].z + temp_z;
         temp_pt.x = Cu.Pt[0].x + temp_x;
		 temp_pt.y = Cu.Pt[0].y + temp_y;
		 temp_pt.z = Cu.Pt[0].z + temp_z;
		 Cu_Backup.AddPt(temp_pt);

		 bool boundary_label = Cu_Backup.Pt[k].check_out_of_range_3D(SM,SN,SZ);

		 index2[0] = Cu_Backup.Pt[k].x;
	     index2[1] = Cu_Backup.Pt[k].y;
	     index2[2] = Cu_Backup.Pt[k].z;
	     v2 = IM->V1->GetPixel(index2);

		 if( boundary_label )
		 {
		    weight(k) = 0;
		 } 
		 else
		 {
		   if( IM->VBW->GetPixel(index2) == 0 || IM->IL_Tracing->GetPixel(index2) != 0 || IM->IL->GetPixel(index2) != 0 )
		   {
		     weight(k) = 0;
		   }
		   else
		   {
		    int v_angle;
            v_angle = acos(v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])/(v1.GetNorm() * v2.GetNorm()) * 180/pi;
		    if( v_angle > 90 )
			   v_angle = 180 - v_angle;
			weight(k) = exp( (float)(90 - v_angle)/60 );
		   }
		   
		 }
	     k++;
		}
	  } 
	}

	//find maximal point and replace tail
	if( weight.max_value() != 0 )
	{
		int idx = weight.arg_max();
	    Cu.Pt[0] = Cu_Backup.Pt[idx];
		jump_gap = true;
	}
	   
   }
   else if( ht == 1 )
   {
 	index1[0] = Cu.Pt[Cu.NP-1].x;
	index1[1] = Cu.Pt[Cu.NP-1].y;
	index1[2] = Cu.Pt[Cu.NP-1].z;
	v1 = IM->V1->GetPixel(index1);

	x = Cu.Pt[Cu.NP-1].x - Cu.Pt[Cu.NP-2].x;
    y = Cu.Pt[Cu.NP-1].y - Cu.Pt[Cu.NP-2].y;
	z = Cu.Pt[Cu.NP-1].z - Cu.Pt[Cu.NP-2].z;
    float theta = acos(z/sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

	float phi = atan2(y,x);

	int k = 0;

	int step_size = 2*angle/step;

	vnl_vector<float> weight(10000);
	weight.fill(0);

	for( int j = 1; j <= dist; j++ )
	{
      for( int i = -angle; i <= angle; i += step_size)
	  {
	    for( int z = -angle/ratio; z <= angle/ratio; z += (step_size/ratio))
		{
		 temp_x = j * sin(theta + (float)z * pi/180) * cos(phi + (float)i * pi/180);
	     temp_y = j * sin(theta + (float)z * pi/180) * sin(phi + (float)i * pi/180);
	     temp_z = j * cos(theta + (float)z * pi/180);
	     //Cu_Backup.Pt[k].x = Cu.Pt[Cu.NP-1].x + temp_x;
         //Cu_Backup.Pt[k].y = Cu.Pt[Cu.NP-1].y + temp_y;
	     //Cu_Backup.Pt[k].z = Cu.Pt[Cu.NP-1].z + temp_z;
		 temp_pt.x = Cu.Pt[Cu.NP-1].x + temp_x;
		 temp_pt.y = Cu.Pt[Cu.NP-1].y + temp_y;
		 temp_pt.z = Cu.Pt[Cu.NP-1].z + temp_z;
		 Cu_Backup.AddPt(temp_pt);

		 bool boundary_label = Cu_Backup.Pt[k].check_out_of_range_3D(SM,SN,SZ);

		 index2[0] = Cu_Backup.Pt[k].x;
	     index2[1] = Cu_Backup.Pt[k].y;
	     index2[2] = Cu_Backup.Pt[k].z;
	     v2 = IM->V1->GetPixel(index2);

		 if( boundary_label )
		 {
		    weight(k) = 0;
		 } 
		 else
		 {
		   if( IM->VBW->GetPixel(index2) == 0 || IM->IL_Tracing->GetPixel(index2) != 0 || IM->IL->GetPixel(index2) != 0 )
		   {
		     weight(k) = 0;
		   }
		   else
		   {
		    int v_angle;
            v_angle = acos(v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])/(v1.GetNorm() * v2.GetNorm()) * 180/pi;
		    if( v_angle > 90 )
			   v_angle = 180 - v_angle;
			weight(k) = exp( (float)(90 - v_angle)/60 );
		   }
		   
		 }
	     k++;
		}
	  } 
	}

	//find maximal point and replace head
	if( weight.max_value() != 0 )
	{
		int idx = weight.arg_max();
	    Cu.Pt[Cu.NP-1] = Cu_Backup.Pt[idx];
		jump_gap = true;
	}
   }
  return jump_gap;
}

bool SnakeClass::Jump_Over_Crossover_New(int dist, int collision)
{

    bool jump_crossover = false;
	int crossover_dist = dist;
    int angle_th = 120;

	//find nearby tip point of existing snake
	if( collision == 1 || collision == 3 )
	{
	  for( int i = 0; i < SnakeList->NSnakes; i++ )
	  {
		  //if( SnakeList->valid_list(i) == 0 || i == tail_collision_snake_id )
		  if( SnakeList->valid_list(i) == 0 )
			  continue;

		  float dist_head = Cu.GetFirstPt().GetDistTo( SnakeList->Snakes[i].Cu.GetLastPt() );
		  float dist_tail = Cu.GetFirstPt().GetDistTo( SnakeList->Snakes[i].Cu.GetFirstPt() );
		  //tail-head merging
		  if( dist_head < dist_tail && dist_head < crossover_dist )
		  {
		    int pt_idx = SnakeList->Snakes[i].Cu.GetSize() - 1;
		  	 //transfer points
			 Cu_Backup = Cu;
            for( int j = pt_idx; j >= 0; j-- )
	        {
		      Cu.AddTailPt( SnakeList->Snakes[i].Cu.Pt[j] );
	        }
		    if( Cu.check_for_sharp_turn(angle_th) )
			{
			  Cu = Cu_Backup;
			}
			else
			{
			 SnakeList->RemoveSnake(i);
			 jump_crossover = true;
			}
			//return jump_crossover;
		  }
		  //tail-tail merging
		  else if( dist_tail < dist_head && dist_tail < crossover_dist )
		  {
		    int pt_idx = 0;
		  	//transfer points
			 Cu_Backup = Cu;
            for( int j = pt_idx; j < SnakeList->Snakes[i].Cu.GetSize(); j++ )
	        {
		      Cu.AddTailPt( SnakeList->Snakes[i].Cu.Pt[j] );
	        }
		    if( Cu.check_for_sharp_turn(angle_th) )
			{
			  Cu = Cu_Backup;
			}
			else
			{
			 SnakeList->RemoveSnake(i);
			 jump_crossover = true;
			}
			//return jump_crossover;
		  }
	  }
	}

	if( collision == 2 || collision == 3 )
	{
	  for( int i = 0; i < SnakeList->NSnakes; i++ )
	  {
		  //if( SnakeList->valid_list(i) == 0 || i == head_collision_snake_id )
		  if( SnakeList->valid_list(i) == 0 )
			  continue;

		  float dist_head = Cu.GetLastPt().GetDistTo( SnakeList->Snakes[i].Cu.GetLastPt() );
		  float dist_tail = Cu.GetLastPt().GetDistTo( SnakeList->Snakes[i].Cu.GetFirstPt() );
		  //head-head merging
		  if( dist_head < dist_tail && dist_head < crossover_dist)
		  {
		    int pt_idx = SnakeList->Snakes[i].Cu.GetSize() - 1;
		  	 //transfer points
			 Cu_Backup = Cu;
            for( int j = pt_idx; j >= 0; j-- )
	        {
		      Cu.AddPt( SnakeList->Snakes[i].Cu.Pt[j] );
	        }
		    if( Cu.check_for_sharp_turn(angle_th) )
			{
			  Cu = Cu_Backup;
			}
			else
			{
			 SnakeList->RemoveSnake(i);
			 jump_crossover = true;
			}
			//return jump_crossover;
		  }
		  //head-tail merging
		  else if( dist_tail < dist_head && dist_tail < crossover_dist )
		  {
		    int pt_idx = 0;
		  	//transfer points
            Cu_Backup = Cu;
            for( int j = pt_idx; j < SnakeList->Snakes[i].Cu.GetSize(); j++ )
	        {
		      Cu.AddPt( SnakeList->Snakes[i].Cu.Pt[j] );
	        }
		    if( Cu.check_for_sharp_turn(angle_th) )
			{
			  Cu = Cu_Backup;
			}
			else
			{
			 SnakeList->RemoveSnake(i);
			 jump_crossover = true;
			}
			//return jump_crossover;
		  }
	  }
	}
	return jump_crossover;
}

bool SnakeClass::Jump_Over_Crossover(int dist, int angle, int step, int ratio, int collision)
{
	bool jump_crossover = false;

	int SM = IM->SM;
    int SN = IM->SN;
    int SZ = IM->SZ;

	//Cu_Backup.NP = 0;
	Cu_Backup.RemoveAllPts();

	int pi = 3.14159265;
	float x,y,z,temp_x,temp_y,temp_z;
    GradientImageType::IndexType index1,index2;
	GradientPixelType v1,v2;

   if( collision == 1 )
   {
	index1[0] = Cu.Pt[0].x;
	index1[1] = Cu.Pt[0].y;
	index1[2] = Cu.Pt[0].z;
	v1 = IM->V1->GetPixel(index1);

	x = Cu.Pt[0].x - Cu.Pt[1].x;
    y = Cu.Pt[0].y - Cu.Pt[1].y;
	z = Cu.Pt[0].z - Cu.Pt[1].z;
    float theta = acos(z/sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

	float phi = atan2(y,x);

	int k = 0;

	int step_size = 2*angle/step;

	vnl_vector<float> weight(10000);
	weight.fill(0);

	for( int j = dist; j <= dist * 2; j++ )
	{
      for( int i = -angle; i <= angle; i += step_size)
	  {
	    for( int z = -angle/ratio; z <= angle/ratio; z += (step_size/ratio))
		{
		 temp_x = j * sin(theta + (float)z * pi/180) * cos(phi + (float)i * pi/180);
	     temp_y = j * sin(theta + (float)z * pi/180) * sin(phi + (float)i * pi/180);
	     temp_z = j * cos(theta + (float)z * pi/180);
	     //Cu_Backup.Pt[k].x = Cu.Pt[0].x + temp_x;
         //Cu_Backup.Pt[k].y = Cu.Pt[0].y + temp_y;
	     //Cu_Backup.Pt[k].z = Cu.Pt[0].z + temp_z;
		 temp_pt.x = Cu.Pt[0].x + temp_x;
		 temp_pt.y = Cu.Pt[0].y + temp_y;
		 temp_pt.z = Cu.Pt[0].z + temp_z;
		 Cu_Backup.AddPt(temp_pt);

		 bool boundary_label = Cu_Backup.Pt[k].check_out_of_range_3D(SM,SN,SZ);

		 index2[0] = Cu_Backup.Pt[k].x;
	     index2[1] = Cu_Backup.Pt[k].y;
	     index2[2] = Cu_Backup.Pt[k].z;
	     v2 = IM->V1->GetPixel(index2);

			 if( !boundary_label && IM->VBW->GetPixel(index2) != 0 && IM->IL_Tracing->GetPixel(index2) == 0 && IM->IL->GetPixel(index2) == 0 )
			 {
			   int v_angle;
			   v_angle = acos(v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])/(v1.GetNorm() * v2.GetNorm()) * 180/pi;
               Cu.Pt[0] = Cu_Backup.Pt[k];
			   jump_crossover = true;
			   return jump_crossover;
			 }
		 
	     k++;
	     //Cu_Backup.NP++;
		}
	  } 
	}

   }
   else if( collision == 2 )
   {
 	index1[0] = Cu.Pt[Cu.NP-1].x;
	index1[1] = Cu.Pt[Cu.NP-1].y;
	index1[2] = Cu.Pt[Cu.NP-1].z;
	v1 = IM->V1->GetPixel(index1);

	x = Cu.Pt[Cu.NP-1].x - Cu.Pt[Cu.NP-2].x;
    y = Cu.Pt[Cu.NP-1].y - Cu.Pt[Cu.NP-2].y;
	z = Cu.Pt[Cu.NP-1].z - Cu.Pt[Cu.NP-2].z;
    float theta = acos(z/sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

	float phi = atan2(y,x);

	int k = 0;

	int step_size = 2*angle/step;

	vnl_vector<float> weight(10000);
	weight.fill(0);

	for( int j = 1; j <= dist; j++ )
	{
      for( int i = -angle; i <= angle; i += step_size)
	  {
	    for( int z = -angle/ratio; z <= angle/ratio; z += (step_size/ratio))
		{
		 temp_x = j * sin(theta + (float)z * pi/180) * cos(phi + (float)i * pi/180);
	     temp_y = j * sin(theta + (float)z * pi/180) * sin(phi + (float)i * pi/180);
	     temp_z = j * cos(theta + (float)z * pi/180);
	     //Cu_Backup.Pt[k].x = Cu.Pt[Cu.NP-1].x + temp_x;
         //Cu_Backup.Pt[k].y = Cu.Pt[Cu.NP-1].y + temp_y;
	     //Cu_Backup.Pt[k].z = Cu.Pt[Cu.NP-1].z + temp_z;
		 temp_pt.x = Cu.Pt[Cu.NP-1].x + temp_x;
		 temp_pt.y = Cu.Pt[Cu.NP-1].y + temp_y;
		 temp_pt.z = Cu.Pt[Cu.NP-1].z + temp_z;
		 Cu_Backup.AddPt(temp_pt);

		 bool boundary_label = Cu_Backup.Pt[k].check_out_of_range_3D(SM,SN,SZ);

		 index2[0] = Cu_Backup.Pt[k].x;
	     index2[1] = Cu_Backup.Pt[k].y;
	     index2[2] = Cu_Backup.Pt[k].z;
	     v2 = IM->V1->GetPixel(index2);

			 if( !boundary_label && IM->VBW->GetPixel(index2) != 0 && IM->IL_Tracing->GetPixel(index2) == 0 && IM->IL->GetPixel(index2) == 0 )
			 {
			   int v_angle;
			   v_angle = acos(v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])/(v1.GetNorm() * v2.GetNorm()) * 180/pi;
               Cu.Pt[Cu.NP-1] = Cu_Backup.Pt[k];
			   jump_crossover = true;
			   return jump_crossover;
			 }

	     k++;
	     //Cu_Backup.NP++;
		}
	  } 
	}
   }

  return jump_crossover;
}


bool SnakeClass::Check_Validity(float minimum_length, float repeat_ratio, int repeat_dist, int snake_id)
{
   int SM = IM->SM;
   int SN = IM->SN;
   int SZ = IM->SZ;

   bool valid = true;

   //float repeat_ratio = 0.3;
   //int repeat_dist = 3;

   float background_percent = 0.5;

   //Point3D temp_pt;
   
   //check repeat tracing
   int label_tail = 0;
   int label_head = 0;

   int overlap_pt = 0;
   for( int i = 0; i < Cu.GetSize(); i++ )
   {

    LabelImageType::IndexType index;
	index[0] = Cu.Pt[i].x;
	index[1] = Cu.Pt[i].y;
	index[2] = Cu.Pt[i].z;
    bool overlap = false;
    for( int ix = -repeat_dist; ix <= repeat_dist; ix++ )
    {
     for( int iy = -repeat_dist; iy <= repeat_dist; iy++ )
	 {
	   for( int iz = -repeat_dist; iz <= repeat_dist; iz++ )
	   { 
          LabelImageType::IndexType new_index;
	      temp_pt.x = index[0] + ix;
          temp_pt.y = index[1] + iy;
		  temp_pt.z = index[2] + iz;
		  temp_pt.check_out_of_range_3D(SM,SN,SZ);
		  new_index[0] = temp_pt.x;
		  new_index[1] = temp_pt.y;
		  new_index[2] = temp_pt.z;
		  
          if( i == 0 )
		  {
		    int label_temp = IM->IL->GetPixel( new_index );
			if( label_temp != 0 )
		       label_tail = label_temp;
		  }
		  else if( i == Cu.GetSize() - 1 )
		  { 
		    int label_temp = IM->IL->GetPixel( new_index );
			if( label_temp != 0 )
		       label_head = label_temp;
		  }

		  if( IM->IL->GetPixel( new_index ) != 0 && IM->IL->GetPixel( new_index ) != snake_id ) 
            overlap = true;
	   }
	 }
	}
	if( overlap )
     overlap_pt++;

   }


   if( label_tail != 0 && label_head != 0 && label_tail != label_head )
   {
	    valid = true;
		return valid;
   }

  //check self intersect
  if( Cu.NP > 6 )
  {
   for( int i = 5; i < Cu.NP; i++ )
   {
	 if( Cu.Pt[0].GetDistTo(Cu.Pt[i]) < repeat_dist )
	 {
	   std::cout<<"self intersection"<<std::endl;
	   valid = false;
	   return valid;
	 }
   }
   for( int i = Cu.NP-6; i >= 0; i-- )
   {
     if( Cu.Pt[Cu.NP-1].GetDistTo(Cu.Pt[i]) < repeat_dist )
	 {
	   std::cout<<"self intersection"<<std::endl;
	   valid = false;
	   return valid;
	 }
   }
  }

   //check length
   //if( Cu.GetLength() <= minimum_length && label_tail == 0 && label_head == 0)
   if( Cu.GetLength() <= minimum_length )
   {
	   std::cout<<"less than minimum length"<<std::endl;
	   valid = false;
	   return valid;
   }

   //check for repeating and loop
   //if( (float)overlap_pt/(float)Cu.GetSize() > repeat_ratio  || (label_tail == label_head && label_tail != 0) )
   if( (float)overlap_pt/(float)Cu.GetSize() > repeat_ratio )
   {
	  std::cout<<"repeating tracing"<<std::endl;
      valid = false;
	  return valid;     
   }

   /*//check for background leakage
   int leak_pt = 0;
   for( int i = 0; i < Cu.GetSize(); i++ )
   {
	   ImageType::IndexType temp;
	   temp[0] = Cu.Pt[i].x;
	   temp[1] = Cu.Pt[i].y;
	   temp[2] = Cu.Pt[i].z;
	   if( IM->VBW->GetPixel(temp) == 0 )
	   {
	     leak_pt++;
	   }
   }
   if( (float)leak_pt/(float)Cu.GetSize() > background_percent )
   {
	   std::cout<<"leaked tracing"<<std::endl;
	   valid = false;
	   return valid;
   }*/

   return valid;
}

void SnakeClass::SetTracedSnakes(SnakeListClass *S)
{
   SnakeList = S;
}

bool SnakeClass::Check_Head_Collision(ImageType::IndexType index, int collision_dist, int minimum_length, bool automatic_merging, int max_angle, int snake_id)
{
   int angle_th = max_angle;

   bool merging = false;
   if( collision == 2 || collision == 3 )
	   return merging;
  
   int SM = IM->I->GetLargestPossibleRegion().GetSize()[0];
   int SN = IM->I->GetLargestPossibleRegion().GetSize()[1];
   int SZ = IM->I->GetLargestPossibleRegion().GetSize()[2];
 
   bool overlap = false;
   //Point3D temp_pt;
   //Point3D head_pt(index[0],index[1],index[2]);
   head_pt.x = index[0];
   head_pt.y = index[1];
   head_pt.z = index[2];
   float L3 = Cu.GetLength();

   for( int ix = -collision_dist; ix <= collision_dist; ix++ )
   {
     for( int iy = -collision_dist; iy <= collision_dist; iy++ )
	 {
	   for( int iz = -collision_dist; iz <= collision_dist; iz++ )
	   { 
          LabelImageType::IndexType new_index;
	      temp_pt.x = index[0] + ix;
          temp_pt.y = index[1] + iy;
		  temp_pt.z = index[2] + iz;
		  temp_pt.check_out_of_range_3D(SM,SN,SZ);
		  new_index[0] = temp_pt.x;
		  new_index[1] = temp_pt.y;
		  new_index[2] = temp_pt.z;


		 //check if the head is in soma region
	     if( IM->Centroid.NP != 0 )
		 {
		  int id_soma = IM->ISoma->GetPixel(new_index);
		  if( id_soma != 0 && Cu.NP > 3 && Cu.GetLength() > minimum_length )
		  {
		     Cu.AddPt( IM->Centroid.Pt[id_soma-1] );
			 Ru.push_back( 0 );
             //BranchPt.AddPt(Cu.Pt[Cu.NP-1]);

			 std::cout<<"head grows into soma region..."<<std::endl;
			 RootPt.AddPt(Cu.Pt[Cu.NP-1]);

			   if( collision == 1 )
				 collision = 3;
			   else
				 collision = 2;
			   merging = true;
			   return merging;
		  }
		 }


		  int id = IM->IL->GetPixel( new_index );
	      if( id != 0 && id != snake_id )
		  {
		    //std::cout<<"head collision"<<std::endl;
			//int id = IM->IL->GetPixel( new_index );

			//find nearest point at the traced snake
			vnl_vector<float> dist_temp(SnakeList->Snakes[id-1].Cu.GetSize());
		    for( int i = 0; i < SnakeList->Snakes[id-1].Cu.GetSize(); i++ )
			{
			   dist_temp(i) =  head_pt.GetDistTo(SnakeList->Snakes[id-1].Cu.Pt[i]);
			}

			int pt_id = dist_temp.arg_min();


			if( automatic_merging )
		    {
            float L1 = SnakeList->Snakes[id-1].Cu.GetPartLength(pt_id,0);
		    float L2 = SnakeList->Snakes[id-1].Cu.GetPartLength(pt_id,1);
			
			if( pt_id != 0 && pt_id != SnakeList->Snakes[id-1].Cu.GetSize() - 1 )
			{
	        
			  if( L1 > minimum_length && L2 > minimum_length && L3 > minimum_length )
			  {

			     Cu.Pt[Cu.NP-1] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];
			     BranchPt.AddPt(Cu.Pt[Cu.NP-1]);
			     std::cout<<"head collision, branch point detected..."<<std::endl;
				 if( collision == 1 )
				  collision = 3;
			     else
				  collision = 2;
				 head_collision_snake_id = id-1;

			     return merging;
			  }
			  else if( L1 <= minimum_length && L3 > minimum_length )
			  {
			     
				 //PointList3D Cu_Temp;
				 Cu_Backup = Cu;
				 std::vector<float> Ru_Backup = Ru;

				 for( int im = pt_id; im < SnakeList->Snakes[id-1].Cu.GetSize(); im++ )
				 {
					 Cu.AddPt(SnakeList->Snakes[id-1].Cu.Pt[im]);
					 Ru.push_back(SnakeList->Snakes[id-1].Ru[im]);
				 }
				 
				 if( Cu.check_for_sharp_turn(angle_th) )
				 {
				   Cu = Cu_Backup;
				   Ru = Ru_Backup;

				   Cu.Pt[Cu.NP-1] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];
			       BranchPt.AddPt(Cu.Pt[Cu.NP-1]);
			       std::cout<<"head collision, branch point detected..."<<std::endl;
				   if( collision == 1 )
				    collision = 3;
			       else
				    collision = 2;
				   head_collision_snake_id = id-1;

				   return merging;
				 }
                 
				 std::cout<<"head merging 1.................."<<std::endl;
				 SnakeList->RemoveSnake(id-1);
				 merging = true;
				 //std::cout<<"Ru size:"<<Ru.size()<<std::endl;
			     //std::cout<<"Cu size:"<<Cu.NP<<std::endl;
				 return merging; 
			  }
			  else if( L2 <= minimum_length && L3 > minimum_length )
			  {
                 
				 //PointList3D Cu_Temp;
				 Cu_Backup = Cu;
                 std::vector<float> Ru_Backup = Ru;

				 for( int im = pt_id; im >= 0; im-- )
				 {
					 Cu.AddPt(SnakeList->Snakes[id-1].Cu.Pt[im]);
					 Ru.push_back(SnakeList->Snakes[id-1].Ru[im]);
				 }

				 if( Cu.check_for_sharp_turn(angle_th) )
				 {
				   Cu = Cu_Backup;
				   Ru = Ru_Backup;

				   Cu.Pt[Cu.NP-1] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];
			       BranchPt.AddPt(Cu.Pt[Cu.NP-1]);
			       std::cout<<"head collision, branch point detected..."<<std::endl;
				   if( collision == 1 )
				    collision = 3;
			       else
				    collision = 2;
				   head_collision_snake_id = id-1;

				   return merging;
				 }
                 
				 std::cout<<"head merging 2.................."<<std::endl;
				 SnakeList->RemoveSnake(id-1);
				 merging = true;
				 //std::cout<<"Ru size:"<<Ru.size()<<std::endl;
			     //std::cout<<"Cu size:"<<Cu.NP<<std::endl;
				 return merging;
			  }
			}
			else if( pt_id == 0 && L3 > minimum_length )
			{
               
			   //PointList3D Cu_Temp;
			   Cu_Backup = Cu;
			   std::vector<float> Ru_Backup = Ru;

			  // Cu.AddPtList(SnakeList->Snakes[id-1].Cu);
			   for( int im = 0; im < SnakeList->Snakes[id-1].Cu.NP; im++ )
			   {
			     Cu.AddPt(SnakeList->Snakes[id-1].Cu.Pt[im]);
			     Ru.push_back(SnakeList->Snakes[id-1].Ru[im]);
			   }

			    if( Cu.check_for_sharp_turn(angle_th) )
				 {
				   Cu = Cu_Backup;
				   Ru = Ru_Backup;

				   Cu.Pt[Cu.NP-1] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];
			       BranchPt.AddPt(Cu.Pt[Cu.NP-1]);
			       std::cout<<"head collision, branch point detected..."<<std::endl;
				   if( collision == 1 )
				    collision = 3;
			       else
				    collision = 2;
				   head_collision_snake_id = id-1;

				   return merging;
				 }

			   std::cout<<"head merging 3.................."<<std::endl;
			   SnakeList->RemoveSnake(id-1);
			   merging = true;
			   //std::cout<<"Ru size:"<<Ru.size()<<std::endl;
			   //std::cout<<"Cu size:"<<Cu.NP<<std::endl;
			   return merging;
			}
			else if( pt_id == SnakeList->Snakes[id-1].Cu.GetSize() - 1 && L3 > minimum_length )
			{
			   SnakeList->Snakes[id-1].Cu.Flip();
               
			   //PointList3D Cu_Temp;
			   Cu_Backup = Cu;
			   std::vector<float> Ru_Backup = Ru;

			   //Cu.AddPtList(SnakeList->Snakes[id-1].Cu);
			   for( int im = 0; im < SnakeList->Snakes[id-1].Cu.NP; im++ )
			   {
			     Cu.AddPt(SnakeList->Snakes[id-1].Cu.Pt[im]);
			     Ru.push_back(SnakeList->Snakes[id-1].Ru[im]);
			   }

			    if( Cu.check_for_sharp_turn(angle_th) )
				 {
				   Cu = Cu_Backup;
				   Ru = Ru_Backup;

				   SnakeList->Snakes[id-1].Cu.Flip();

				   Cu.Pt[Cu.NP-1] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];
			       BranchPt.AddPt(Cu.Pt[Cu.NP-1]);
			       std::cout<<"head collision, branch point detected..."<<std::endl;
				   if( collision == 1 )
				    collision = 3;
			       else
				    collision = 2;
				   head_collision_snake_id = id-1;

				   return merging;
				 }

			   std::cout<<"head merging 4.................."<<std::endl;
			   SnakeList->RemoveSnake(id-1);
			   merging = true;
			   //std::cout<<"Ru size:"<<Ru.size()<<std::endl;
			   //std::cout<<"Cu size:"<<Cu.NP<<std::endl;
			   return merging;
			}
			}
			else
			{
              Cu.Pt[Cu.NP-1] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];
			  BranchPt.AddPt(Cu.Pt[Cu.NP-1]);
              //SnakeList->Snakes[id-1].BranchPt.AddPt(Cu.Pt[Cu.NP-1]);
			  merging = false;
			  std::cout<<"head collision, branch point detected..."<<std::endl;
			  if( collision == 1 )
				collision = 3;
			  else
				collision = 2;
			  head_collision_snake_id = id-1;
			  return merging;
            }
		  }
	   }
	 }
   }
   return merging;
}

bool SnakeClass::Check_Tail_Collision(ImageType::IndexType index, int collision_dist, int minimum_length, bool automatic_merging, int max_angle, int snake_id)
{
   int angle_th = max_angle;
   bool merging = false;

   if( collision == 1 || collision == 3)
	   return merging;

   int SM = IM->I->GetLargestPossibleRegion().GetSize()[0];
   int SN = IM->I->GetLargestPossibleRegion().GetSize()[1];
   int SZ = IM->I->GetLargestPossibleRegion().GetSize()[2];
 
   bool overlap = false;
   //Point3D temp_pt;
   //Point3D tail_pt(index[0],index[1],index[2]);
   tail_pt.x = index[0];
   tail_pt.y = index[1];
   tail_pt.z = index[2];

   float L3 = Cu.GetLength();

   for( int ix = -collision_dist; ix <= collision_dist; ix++ )
   {
     for( int iy = -collision_dist; iy <=collision_dist; iy++ )
	 {
	   for( int iz = -collision_dist; iz <=collision_dist; iz++ )
	   { 
          LabelImageType::IndexType new_index;

	      temp_pt.x = index[0] + ix;
          temp_pt.y = index[1] + iy;
		  temp_pt.z = index[2] + iz;
		  temp_pt.check_out_of_range_3D(SM,SN,SZ);
		  new_index[0] = temp_pt.x;
		  new_index[1] = temp_pt.y;
		  new_index[2] = temp_pt.z;

		 //check if the tail is in soma region
	     if( IM->Centroid.NP != 0 )
		 {
		  int id_soma = IM->ISoma->GetPixel(new_index);
		  if( id_soma != 0 && Cu.NP > 3 && Cu.GetLength() > minimum_length)
		  {
		     Cu.AddTailPt( IM->Centroid.Pt[id_soma-1] );
			 Ru.push_back( 0 );
			 //BranchPt.AddPt(Cu.Pt[0]);

			 std::cout<<"tail grows into soma region..."<<std::endl;
             RootPt.AddPt(Cu.Pt[0]);

			   if( collision == 2 )
				 collision = 3;
			   else
				 collision = 1;
			   merging = true;
			   return merging;
		  }
		 }

		  int id = IM->IL->GetPixel( new_index );
	      if( id != 0 && id != snake_id )
		  {
			//int id = IM->IL->GetPixel( new_index );
			//find nearest point at the traced snake
			vnl_vector<float> dist_temp(SnakeList->Snakes[id-1].Cu.GetSize());
            //check if snake is removed
			if( SnakeList->valid_list(id-1) == 0 )
                 continue;

			overlap = true;
		    for( int i = 0; i < SnakeList->Snakes[id-1].Cu.GetSize(); i++ )
			{
			   dist_temp(i) =  tail_pt.GetDistTo(SnakeList->Snakes[id-1].Cu.Pt[i]);
			}

			int pt_id = dist_temp.arg_min();

		    if( automatic_merging )
			{
			float L1 = SnakeList->Snakes[id-1].Cu.GetPartLength(pt_id,0);
		    float L2 = SnakeList->Snakes[id-1].Cu.GetPartLength(pt_id,1);

			if( pt_id != 0 && pt_id != SnakeList->Snakes[id-1].Cu.GetSize() - 1 )
			{
	        
			  if( L1 > minimum_length && L2 > minimum_length && L3 > minimum_length )
			  {
	
			     Cu.Pt[0] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];
			     BranchPt.AddPt(Cu.Pt[0]);
			     std::cout<<"tail collision, branch point detected..."<<std::endl;
				 if( collision == 2 )
				  collision = 3;
			     else
				  collision = 1;
				 tail_collision_snake_id = id-1;

			     return merging;
			  }
			  else if( L1 <= minimum_length && L2 > minimum_length && L3 > minimum_length )
			  {

				 SnakeList->Snakes[id-1].Cu.Flip();

                 //PointList3D Cu_Temp;
				 Cu_Backup = Cu;
				 std::vector<float> Ru_Backup = Ru;

				 for( int im = SnakeList->Snakes[id-1].Cu.GetSize() - pt_id; im >= 0; im-- )
				 {
					 Cu.AddTailPt(SnakeList->Snakes[id-1].Cu.Pt[im]);
					 std::vector<float>::iterator it;
                     it = Ru.begin();
                     it = Ru.insert ( it , SnakeList->Snakes[id-1].Ru[im] );
				 }

			    if( Cu.check_for_sharp_turn(angle_th) )
				 {
				   Cu = Cu_Backup;
				   Ru = Ru_Backup;

				   SnakeList->Snakes[id-1].Cu.Flip();

				   Cu.Pt[0] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];
			       BranchPt.AddPt(Cu.Pt[0]);
			       std::cout<<"tail collision, branch point detected..."<<std::endl;
				   if( collision == 2 )
				    collision = 3;
			       else
				    collision = 1;
				   tail_collision_snake_id = id-1;

				   return merging;
				 }

				 std::cout<<"tail merging 1.................."<<std::endl;
				 SnakeList->RemoveSnake(id-1);
				 merging = true;
				 return merging; 
			  }
			  else if( L1 > minimum_length && L2 <= minimum_length && L3 > minimum_length )
			  {
                 //PointList3D Cu_Temp;

				 Cu_Backup = Cu;
                 std::vector<float> Ru_Backup = Ru;

				 for( int im = pt_id; im >= 0; im-- )
				 {
					 Cu.AddTailPt(SnakeList->Snakes[id-1].Cu.Pt[im]);
					 std::vector<float>::iterator it;
                     it = Ru.begin();
                     it = Ru.insert ( it , SnakeList->Snakes[id-1].Ru[im] );
				 }

			    if( Cu.check_for_sharp_turn(angle_th) )
				 {
				   Cu = Cu_Backup;
				   Ru = Ru_Backup;

				   Cu.Pt[0] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];

				   BranchPt.RemoveAllPts();
			       BranchPt.AddPt(Cu.Pt[0]);
			       std::cout<<"tail collision, branch point detected..."<<std::endl;
				   if( collision == 2 )
				    collision = 3;
			       else
				    collision = 1;
				   tail_collision_snake_id = id-1;

				   return merging;
				 }
   
				 std::cout<<"tail merging 2.................."<<std::endl;
				 SnakeList->RemoveSnake(id-1);
				 merging = true;
				 return merging;
			  }
			}
			else if( pt_id == 0 && L3 > minimum_length )
			{
			   SnakeList->Snakes[id-1].Cu.Flip();

			   //PointList3D Cu_Temp;

			   Cu_Backup = Cu;
			   std::vector<float> Ru_Backup = Ru;

			   //Cu.AddTailPtList(SnakeList->Snakes[id-1].Cu);
               for( int im = SnakeList->Snakes[id-1].Cu.GetSize()-1 ; im >=0; im-- )
			   {
			      Cu.AddTailPt(SnakeList->Snakes[id-1].Cu.Pt[im]);
				  std::vector<float>::iterator it;
                  it = Ru.begin();
                  it = Ru.insert ( it , SnakeList->Snakes[id-1].Ru[im] );
			   }

			    if( Cu.check_for_sharp_turn(angle_th) )
				{
				   Cu = Cu_Backup;
				   Ru = Ru_Backup;

				   SnakeList->Snakes[id-1].Cu.Flip();

				   Cu.Pt[0] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];

			       BranchPt.AddPt(Cu.Pt[0]);
			       std::cout<<"tail collision, branch point detected..."<<std::endl;
				   if( collision == 2 )
				    collision = 3;
			       else
				    collision = 1;
				   tail_collision_snake_id = id-1;

				   return merging;
				}

			   std::cout<<"tail merging 3.................."<<std::endl;
			   SnakeList->RemoveSnake(id-1);
			   merging = true;
			   return merging;
			}
			else if( pt_id == SnakeList->Snakes[id-1].Cu.GetSize() - 1 && L3 > minimum_length )
			{
	           //PointList3D Cu_Temp;
			   Cu_Backup = Cu;
			   std::vector<float> Ru_Backup = Ru;

			   //Cu.AddTailPtList(SnakeList->Snakes[id-1].Cu);
               for( int im = pt_id; im >=0; im-- )
			   {
			      Cu.AddTailPt(SnakeList->Snakes[id-1].Cu.Pt[im]);
				  std::vector<float>::iterator it;
                  it = Ru.begin();
                  it = Ru.insert ( it , SnakeList->Snakes[id-1].Ru[im] );
			   }

			    if( Cu.check_for_sharp_turn(angle_th) )
				{
				   Cu = Cu_Backup;
				   Ru = Ru_Backup;

				   Cu.Pt[0] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];

			       BranchPt.AddPt(Cu.Pt[0]);
			       std::cout<<"tail collision, branch point detected..."<<std::endl;
				   if( collision == 2 )
				    collision = 3;
			       else
				    collision = 1;
				   tail_collision_snake_id = id-1;

				   return merging;
				}

			   std::cout<<"tail merging 4.................."<<std::endl;
			   SnakeList->RemoveSnake(id-1);
			   merging = true;
			   return merging;
			}
			}
			else
			{
              Cu.Pt[0] = SnakeList->Snakes[id-1].Cu.Pt[pt_id];
			  BranchPt.AddPt(Cu.Pt[0]);
              // SnakeList->Snakes[id-1].BranchPt.AddPt(Cu.Pt[0]);
			  merging = false;
			  std::cout<<"tail collision, branch point detected..."<<std::endl;
			   if( collision == 2 )
				 collision = 3;
			   else
				 collision = 1;
			   tail_collision_snake_id = id-1;
			   return merging;
			}
		  }
	   }
	 }
   }
   return merging;
}

bool SnakeClass::Check_Head_Leakage(ImageType::IndexType head_index)
{
	 bool leakage = true;
	 for( int ix = -1; ix <= 1; ix++ )
	 {
      for( int iy = -1; iy <= 1; iy++ )
	  {
		for( int iz = -1; iz <= 1; iz++ )
		{
			ImageType::IndexType idx_temp;
			temp_pt.x = head_index[0] + ix;
			temp_pt.y = head_index[1] + iy;
			temp_pt.z = head_index[2] + iz;
			temp_pt.check_out_of_range_3D(IM->SM,IM->SN,IM->SZ);
			idx_temp[0] = temp_pt.x;
            idx_temp[1] = temp_pt.y;
			idx_temp[2] = temp_pt.z;

			int end = Cu.NP - 1;

			int d_product = (Cu.Pt[end-1].x - Cu.Pt[end-2].x) * (idx_temp[0]- Cu.Pt[end-1].x) + 
				            (Cu.Pt[end-1].y - Cu.Pt[end-2].y) * (idx_temp[1]- Cu.Pt[end-1].y) +
							(Cu.Pt[end-1].z - Cu.Pt[end-2].z) * (idx_temp[2]- Cu.Pt[end-1].z);

			if( IM->VBW->GetPixel(idx_temp) != 0 && IM->IL_Tracing->GetPixel(idx_temp) == 0 && d_product > 0 )
			{  
                Cu.Pt[Cu.NP-1].x = idx_temp[0];
                Cu.Pt[Cu.NP-1].y = idx_temp[1];
				Cu.Pt[Cu.NP-1].z = idx_temp[2];
				leakage = false;
				return leakage;
			}
		}
		}
	  }

	 return leakage;
}

bool SnakeClass::Check_Tail_Leakage(ImageType::IndexType tail_index)
{
	 bool leakage = true;
	 for( int ix = -1; ix <= 1; ix++ )
	 {
      for( int iy = -1; iy <= 1; iy++ )
	  {
		for( int iz = -1; iz <= 1; iz++ )
		{
			ImageType::IndexType idx_temp;
			temp_pt.x = tail_index[0] + ix;
			temp_pt.y = tail_index[1] + iy;
			temp_pt.z = tail_index[2] + iz;
			temp_pt.check_out_of_range_3D(IM->SM,IM->SN,IM->SZ);
			idx_temp[0] = temp_pt.x;
            idx_temp[1] = temp_pt.y;
			idx_temp[2] = temp_pt.z;

		    int d_product = (Cu.Pt[1].x - Cu.Pt[2].x) * (idx_temp[0]- Cu.Pt[1].x) + 
				            (Cu.Pt[1].y - Cu.Pt[2].y) * (idx_temp[1]- Cu.Pt[1].y) +
							(Cu.Pt[1].z - Cu.Pt[2].z) * (idx_temp[2]- Cu.Pt[1].z);

			if( IM->VBW->GetPixel(idx_temp) != 0 && IM->IL_Tracing->GetPixel(idx_temp) == 0 && d_product > 0 )
			{  
                Cu.Pt[0].x = idx_temp[0];
                Cu.Pt[0].y = idx_temp[1];
				Cu.Pt[0].z = idx_temp[2];
				leakage = false;
				return leakage;
			}
		}
	  }
	 }

	 return leakage;
}

void SnakeClass::Estimate_Radius()
{ 

 /*if( Cu.NP < 3 )
    return;

 int pi = 3.1415926;
 int N = 4;
 int r_min = 1;
 int r_max = 3;
 float step = 0.5;
 vnl_vector<float> radius( (int)(r_max-r_min)/step + 1 );
 for( int i = 0; i < radius.size(); i++ )
 {
   radius(i) = r_min + step * i;
 }


 Ru.set_size(Cu.NP);
 Ru.fill(0);
 Vector3D v1,v2,v3,vtemp,vtemp1;
 Point3D temp;

 typedef itk::VectorLinearInterpolateImageFunction< 
                       GradientImageType, float >  GradientInterpolatorType;
 GradientInterpolatorType::Pointer interpolator = GradientInterpolatorType::New();
 interpolator->SetInputImage(IM->IG);

 GradientImageType::IndexType index;

 for( int i = 1; i < Cu.NP-1; i++ )
 {
   vnl_vector<float> score(radius.size());
   score.fill(0);

   //if( i == 0 )
   //{
   //  v1.x = Cu.Pt[1].x - Cu.Pt[0].x;
   //  v1.y = Cu.Pt[1].y - Cu.Pt[0].y;
   //	 v1.z = Cu.Pt[1].z - Cu.Pt[0].z;
   // v1.ConvertUnit();
   //}

   v1.x = Cu.Pt[i].x - Cu.Pt[i-1].x;
   v1.y = Cu.Pt[i].y - Cu.Pt[i-1].y;
   v1.z = Cu.Pt[i].z - Cu.Pt[i-1].z;
   v1.ConvertUnit();
  
   v2.x = -v1.z;
   v2.y = 0;
   v2.z = v1.x;
   v2.ConvertUnit();
   v3.x = 1;
   v3.y = -(pow(v1.x,2) + pow(v1.z,2))/(v1.x*v1.y + std::numeric_limits<double>::epsilon());
   v3.z = v1.z/(v1.x + std::numeric_limits<double>::epsilon());
   v3.ConvertUnit();

   //std::cout<<"v1 * v2:"<<v1.GetDProduct(v2)<<std::endl;
   //std::cout<<"v1 * v3:"<<v1.GetDProduct(v3)<<std::endl;
   //std::cout<<"v2 * v3:"<<v2.GetDProduct(v3)<<std::endl;
   for( int j = 0; j < radius.size(); j++ )
   {
     for( int z = 0; z < N; z++ )
	 {
	    float theta = (2 * pi * z)/N;
	    vtemp.x = v2.x * cos(theta) + v3.x * sin(theta);
	    vtemp.y = v2.y * cos(theta) + v3.y * sin(theta);
		vtemp.z = v2.z * cos(theta) + v3.z * sin(theta);
		vtemp.ConvertUnit();

		//std::cout<<"vtemp.GetDProduct(v1):"<<vtemp.GetDProduct(v1)<<std::endl;

		temp = Cu.Pt[i] + vtemp * radius(j);
		//std::cout<<"temp:"<<v2.x * cos(theta) + v3.x * sin(theta)<<","<<temp.x<<","<<temp.y<<","<<temp.z<<std::endl;
		if( isinf( temp.x ) || isinf( temp.y ) || isinf( temp.z ) )
		{
		  continue;
		}
		
		if( temp.check_out_of_range_3D(IM->SM,IM->SN,IM->SZ) )
			continue;

		index[0] = temp.x;
		index[1] = temp.y;
		index[2] = temp.z;
   
		GradientPixelType gvf = interpolator->EvaluateAtIndex(index);
		vtemp1.x = gvf[0];
		vtemp1.y = gvf[1];
		vtemp1.z = gvf[2];
		vtemp1.ConvertUnit();

		float Cz = -1 * ( vtemp.x * vtemp1.x + vtemp.y * vtemp1.y + vtemp.z * vtemp1.z );
		if( Cz < 0 )
			Cz = 0;
		
		float Rz = sqrt(pow(gvf[0],2) + pow(gvf[1],2) + pow(gvf[2],2)) * pow(Cz,2);

		//std::cout<<"Cz Mag:"<<Cz<<",  "<<sqrt(pow(gvf[0],2) + pow(gvf[1],2) + pow(gvf[2],2))<<std::endl;
		score( j ) += Rz;
	 }
   }

   //for( int k = 0; k < score.size(); k++ )
   //{
	//   std::cout<<score(k)<<",";
   //}
   //std::cout<<std::endl;
   int max_id = score.arg_max();
   Ru(i) = radius(max_id);
   //std::cout<<"Ru(i):"<<Ru(i)<<std::endl;
 }
 
 Ru(0) = Ru(1);
 Ru(Ru.size()-1) = Ru(Ru.size()-2);

 
 //radius smoothing
 int cc = 1;
 float R_sd = 0;
 float R_mean = Ru.mean();
 Ru = Ru  * sqrt((double)2);
 for( int i = 0; i < Ru.size(); i++)
 {
   R_sd += pow(Ru(i)-R_mean,2); 
 }
 R_sd = sqrt(R_sd/Ru.size());

 for( int i = 0; i < Ru.size(); i++)
 {
   if( i == 0 )
   {
      Ru(i) = (Ru(i) + Ru(i+1))/2;
   } 
   else if( i == Ru.size() - 1 )
   {
      Ru(i) = (Ru(i) + Ru(i-1))/2;
   }
   else
   {
     if( Ru(i) > R_mean + cc * R_sd || Ru(i) < R_mean - cc * R_sd)
     {
       Ru(i) = (Ru(i-1) + Ru(i+1))/2;
     }
	 else
	 {
	   Ru(i) = (Ru(i-1) + Ru(i) + Ru(i+1))/3;
	 }
   }

 } */
} 

bool SnakeClass::Compute_Seed_Force(int head_tail, int distance)
{
	bool seed_force = false;
	Vector3D SForce;

   if( head_tail == 0 )
   {
    vnl_vector<float> dist_temp(IM->SeedPt.NP);
	dist_temp.fill(10000);

	SForce.x = Cu.Pt[0].x - Cu.Pt[2].x;
    SForce.y = Cu.Pt[0].x - Cu.Pt[2].y;
	SForce.z = Cu.Pt[0].x - Cu.Pt[2].z;

    for( int i = 0; i < IM->SeedPt.NP; i++ )
    {
	   if( IM->visit_label(i) == 1 )
		   continue;
	   
	   //if( Cu.Pt[0].GetDistTo(IM->SeedPt.Pt[i]) <= distance )
	   //{
	     if( ((IM->SeedPt.Pt[i].x - Cu.Pt[0].x) * SForce.x + (IM->SeedPt.Pt[i].y - Cu.Pt[0].y) * SForce.y + (IM->SeedPt.Pt[i].z - Cu.Pt[0].z) * SForce.z) > 0 )
		 {
		    //seed_force = true;
			//Cu.Pt[0] = IM->SeedPt.Pt[i];
			//break;
           dist_temp(i) = Cu.Pt[0].GetDistTo(IM->SeedPt.Pt[i]);
		 }
	   //}
    }
	
	if( dist_temp.min_value() <= distance )
	{
		seed_force = true;
	    Cu.Pt[0] = IM->SeedPt.Pt[dist_temp.arg_min()];
	}
   }
   else
   {
    vnl_vector<float> dist_temp(IM->SeedPt.NP);
	dist_temp.fill(10000);

   SForce.x = Cu.Pt[Cu.NP-1].x - Cu.Pt[Cu.NP-3].x;
    SForce.y = Cu.Pt[Cu.NP-1].x - Cu.Pt[Cu.NP-3].y;
	SForce.z = Cu.Pt[Cu.NP-1].x - Cu.Pt[Cu.NP-3].z;

    for( int i = 0; i < IM->SeedPt.NP; i++ )
    {
	   if( IM->visit_label(i) == 1 )
		   continue;
	   
	   //if( Cu.Pt[Cu.NP-1].GetDistTo(IM->SeedPt.Pt[i]) <= distance )
	   //{
	     if( ((IM->SeedPt.Pt[i].x - Cu.Pt[Cu.NP-1].x) * SForce.x + (IM->SeedPt.Pt[i].y - Cu.Pt[Cu.NP-1].y) * SForce.y + (IM->SeedPt.Pt[i].z - Cu.Pt[Cu.NP-1].z) * SForce.z) > 0 )
		 {
		    //seed_force = true;
			//Cu.Pt[Cu.NP-1] = IM->SeedPt.Pt[i];
			//break;
			 dist_temp(i) = Cu.Pt[Cu.NP-1].GetDistTo(IM->SeedPt.Pt[i]);
		 }
	   //}
    }
	if( dist_temp.min_value() <= distance)
	{
	  seed_force = true;
      Cu.Pt[Cu.NP-1] = IM->SeedPt.Pt[dist_temp.arg_min()];
	}
   }

   return seed_force;
}

void SnakeClass::OpenSnakeStretch_4D(float alpha, int ITER, int pt_distance, float beta, float kappa, float gamma, 
                                float stretchingRatio, int collision_dist, int minimum_length, 
								bool automatic_merging, int max_angle, bool freeze_body, int s_force, int snake_id, int tracing_model)
{

}


void SnakeClass::OpenSnakeStretch(float alpha, int ITER, int pt_distance, float beta, float kappa, float gamma, 
                                float stretchingRatio, int collision_dist, int N_Active, int minimum_length, 
								bool automatic_merging, int max_angle, bool freeze_body, int s_force, int snake_id)
{

   int hit_boundary_dist = 0;
   int crossover_dist = 10;

   //when collision happened for both tail and head, stop stretching
   if( collision == 3)
	   return;

   typedef itk::VectorLinearInterpolateImageFunction< 
                       GradientImageType, float >  GradientInterpolatorType;
   GradientInterpolatorType::Pointer interpolator = GradientInterpolatorType::New();
   interpolator->SetInputImage(IM->IGVF);
   GradientInterpolatorType::Pointer interpolator_V1 = GradientInterpolatorType::New();
   interpolator_V1->SetInputImage(IM->V1);

   int SM = IM->SM;
   int SN = IM->SN;
   int SZ = IM->SZ;
  
   //int N_Active = 3;

   int N = Cu.GetSize();
   vnl_matrix<float> A = makeOpenA(alpha, beta, N);
   vnl_matrix<float> I(N,N);
   I.set_identity();
   vnl_matrix<float> invAI = vnl_matrix_inverse<float>( A + I * gamma);

   for( int i = 0; i < ITER; i++ )
   {
     vnl_vector<float> vfx(N);
     vnl_vector<float> vfy(N);
	 vnl_vector<float> vfz(N);
     vnl_vector<float> evfx(N);
     vnl_vector<float> evfy(N);
	 vnl_vector<float> evfz(N);
     vnl_vector<float> x(N);
     vnl_vector<float> y(N);
	 vnl_vector<float> z(N);
	 vfx.fill(0);
	 vfy.fill(0);
	 vfz.fill(0);
	 evfx.fill(0);
	 evfy.fill(0);
	 evfz.fill(0);

	 GradientImageType::IndexType head_index; 
	 GradientImageType::IndexType tail_index; 

	 for( int j = 0; j < Cu.GetSize(); j++ )
	 {

		x(j) = Cu.GetPt(j).x;
		y(j) = Cu.GetPt(j).y;
        z(j) = Cu.GetPt(j).z;

        GradientImageType::IndexType index; 
        index[0] = (x(j));
		index[1] = (y(j));
		index[2] = (z(j));

	   if( j == 0 )
	   {
	     tail_index[0] = Cu.GetPt(j).x;
		 tail_index[1] = Cu.GetPt(j).y;
		 tail_index[2] = Cu.GetPt(j).z;

	  	 GradientPixelType v1 = interpolator_V1->EvaluateAtIndex(tail_index);
		 evfx(j) = v1[0];
         evfy(j) = v1[1];
		 evfz(j) = v1[2];

		 /*//project field flow to the normal plane
		 GradientPixelType gvf = interpolator->EvaluateAtIndex(index);
         float dot_product = gvf[0] * evfx(j) + gvf[1] * evfy(j) + gvf[2] * evfz(j);
         vfx(j) = gvf[0] - dot_product * evfx(j);
         vfy(j) = gvf[1] - dot_product * evfy(j);;
	 	 vfz(j) = gvf[2] - dot_product * evfz(j);;
         continue; */
	   }
	   if( j == Cu.GetSize() - 1 )
	   {
	   	 head_index[0] = Cu.GetPt(j).x;
		 head_index[1] = Cu.GetPt(j).y;
		 head_index[2] = Cu.GetPt(j).z;

		 GradientPixelType v1 = interpolator_V1->EvaluateAtIndex(head_index);
		 evfx(j) = v1[0];
         evfy(j) = v1[1];
		 evfz(j) = v1[2];

		 /*//project field flow to the normal plane
		 GradientPixelType gvf = interpolator->EvaluateAtIndex(index);
         float dot_product = gvf[0] * evfx(j) + gvf[1] * evfy(j) + gvf[2] * evfz(j);
         vfx(j) = gvf[0] - dot_product * evfx(j);
         vfy(j) = gvf[1] - dot_product * evfy(j);;
	 	 vfz(j) = gvf[2] - dot_product * evfz(j);;
         continue;*/
	   }


	   if( freeze_body )
	   {
	   	if( j >= Cu.GetSize() - N_Active || j <= N_Active ) 
	    {
		 GradientPixelType gvf = interpolator->EvaluateAtIndex(index);
         vfx(j) = gvf[0];
         vfy(j) = gvf[1];
	 	 vfz(j) = gvf[2];
	    }
	   }
	   else
	   {
	     GradientPixelType gvf = interpolator->EvaluateAtIndex(index);
         vfx(j) = gvf[0];
         vfy(j) = gvf[1];
	 	 vfz(j) = gvf[2];
	   }

	 }

	 //freeze tail or head part when collision happens
	 if( collision == 1 )
  	 {
		for( int k = 0; k < 1; k++ )
		{
		   vfx(k) = 0;
		   vfy(k) = 0;
		   vfz(k) = 0;
		}
	 }
	 else if( collision == 2 )
	 {
	 	for( int k = N-1; k < N; k++ )
		{
		   vfx(k) = 0;
		   vfy(k) = 0;
		   vfz(k) = 0;
		}
	 }

	 Vector3D tailForce( (x(0) - x(2)), (y(0) - y(2)), (z(0) - z(2)) );
	 tailForce.ConvertUnit();
     Vector3D tailForce1( evfx(0), evfy(0), evfz(0) );
     
	 if( tailForce.GetDProduct( tailForce1 ) < 0 )
        tailForce1 = tailForce1 * -1;

	Vector3D tForce;
    if( s_force == 0 )
	   tForce = tailForce1;
	else if( s_force == 1 )
       tForce = tailForce;
	else if( s_force == 2 )
	{
       tForce = (tailForce + tailForce1)/2;
	   tForce.ConvertUnit();
	}

	 int end = x.size() - 1;
     Vector3D headForce( (x(end) - x(end-2)), (y(end) - y(end-2)), (z(end) - z(end-2)) );
	 headForce.ConvertUnit();
     Vector3D headForce1( evfx(end), evfy(end), evfz(end) );
     
	 if( headForce.GetDProduct( headForce1 ) < 0 )
        headForce1 = headForce1 * -1;


	Vector3D hForce;
	if( s_force == 0 )
	   hForce = headForce1;
	else if( s_force == 1 )
       hForce = headForce;
	else if( s_force == 2 )
	{
       hForce = (headForce + headForce1)/2;
	   hForce.ConvertUnit();
	}

	 //check for tail leakage and self-intersection
	 if( IM->IL_Tracing->GetPixel(tail_index) == 1 )
	 {
	     tForce.x = 0;
		 tForce.y = 0;
		 tForce.z = 0;
	 }
	 //if( IM->VBW->GetPixel(tail_index) != snake_id )
	 if( IM->VBW->GetPixel(tail_index) == 0 )
	 {	    
		 tForce.x = 0;
		 tForce.y = 0;
		 tForce.z = 0;
	 }

     //check for tail collision
	 int draw_force_tail = 1;
	
	 bool tail_merging = Check_Tail_Collision(tail_index,collision_dist,minimum_length,automatic_merging,max_angle,snake_id);

	 if( tail_merging )
		break;

	 if( collision == 1 || collision == 3 )
	 {
		x(0) = Cu.Pt[0].x;
		y(0) = Cu.Pt[0].y;
		z(0) = Cu.Pt[0].z;
		tail_index[0] = x(0);
		tail_index[1] = y(0);
		tail_index[2] = z(0);
	    draw_force_tail = 0;
	 }

	 //check for boundary condition
	 int boundary_tail = 0;
	 if( ceil((float)tail_index[0]) >= SM-hit_boundary_dist || ceil((float)tail_index[1]) >= SN-hit_boundary_dist || ceil((float)tail_index[2]) >= SZ-hit_boundary_dist
		 || ceil((float)tail_index[0]) < hit_boundary_dist || ceil((float)tail_index[1]) < hit_boundary_dist || ceil((float)tail_index[2]) < hit_boundary_dist )
	 {
	    draw_force_tail = 0;
	    boundary_tail = 1;
	 }

     //check for head leakage and self_intersection
	 if( IM->IL_Tracing->GetPixel(head_index) == 1 )
	 {
	 	hForce.x = 0;
		hForce.y = 0;
		hForce.z = 0;
	 }

	 //if( IM->VBW->GetPixel(head_index) != snake_id )
	 if( IM->VBW->GetPixel(head_index) == 0 )
	 {
	   //bool head_leakage = Check_Head_Leakage(head_index);
	    hForce.x = 0;
		hForce.y = 0;
		hForce.z = 0;
	 }

	 //check for head collision
	 int draw_force_head = 1;

     bool head_merging = Check_Head_Collision(head_index,collision_dist,minimum_length,automatic_merging,max_angle,snake_id);

	 if( head_merging )
		break;

	 if( collision == 2 || collision == 3 )
	 {
	    x(end) = Cu.Pt[end].x;
		y(end) = Cu.Pt[end].y;
		z(end) = Cu.Pt[end].z;
	    head_index[0] = x(end);
		head_index[1] = y(end);
		head_index[2] = z(end);
	    draw_force_head = 0;
	 }

	 //check for boundary condition
	 int boundary_head = 0;
	 if( ceil((float)head_index[0]) >= SM-hit_boundary_dist || ceil((float)head_index[1]) >= SN-hit_boundary_dist || ceil((float)head_index[2]) >= SZ-hit_boundary_dist
		 || ceil((float)head_index[0]) < hit_boundary_dist || ceil((float)head_index[1]) < hit_boundary_dist || ceil((float)head_index[2]) < hit_boundary_dist )
	 {
	    draw_force_head = 0;
	    boundary_head = 1;
	 }

	//check for crossover
   /*if( jumping_crossover )
   {
	  bool jump_crossover_label = false;
	  //bool jump_crossover_label = Jump_Over_Crossover(crossover_dist, 45, jumping_step, jumping_ratio, collision);
     if( !jump_crossover_label )
	    jump_crossover_label = Jump_Over_Crossover_New(crossover_dist, collision);
    
     if( jump_crossover_label )
	 {
		 std::cout<<"crossover detected"<<std::endl;
		 //BranchPt.NP = 0;
		 BranchPt.RemoveAllPts();
		 collision = 0;
		 return;
	 }

	//bool jump_crossover_label = Jump_Over_Crossover(5, 30, 5, 2, collision);
	if( collision == 1 && jump_crossover_label )
	{
		  std::cout<<"tail jumping over crossover"<<std::endl;
          x(0) = Cu.Pt[0].x;
		  y(0) = Cu.Pt[0].y;
		  z(0) = Cu.Pt[0].z;
		  tail_index[0] = x(0);
		  tail_index[1] = y(0);
		  tail_index[2] = z(0);
		  
		  GradientPixelType gvf = interpolator->EvaluateAtIndex(tail_index);
          vfx(0) = gvf[0];
          vfy(0) = gvf[1];
		  vfz(0) = gvf[2];

		  BranchPt.NP = 0;
		  collision = 0;
	}
	else if( collision == 2 && jump_crossover_label )
	{
		  std::cout<<"head jumping over crossover"<<std::endl;
          x(end) = Cu.Pt[end].x;
		  y(end) = Cu.Pt[end].y;
		  z(end) = Cu.Pt[end].z;
		  head_index[0] = x(end);
		  head_index[1] = y(end);
		  head_index[2] = z(end);
		  
		  GradientPixelType gvf = interpolator->EvaluateAtIndex(head_index);
          vfx(end) = gvf[0];
          vfy(end) = gvf[1];
		  vfz(end) = gvf[2];

		  BranchPt.NP = 0;
		  collision = 0;
	}
   }*/

    vfx(0) = (vfx(0) + stretchingRatio * tForce.x) * draw_force_tail;
    vfy(0) = (vfy(0) + stretchingRatio * tForce.y) * draw_force_tail;
    vfz(0) = (vfz(0) + stretchingRatio * tForce.z) * draw_force_tail;  
   
    vfx(end) = (vfx(end) + stretchingRatio * hForce.x) * draw_force_head;
    vfy(end) = (vfy(end) + stretchingRatio * hForce.y) * draw_force_head;
    vfz(end) = (vfz(end) + stretchingRatio * hForce.z) * draw_force_head;


    x = (invAI * ( x * gamma +  vfx * kappa));
    y = (invAI * ( y * gamma +  vfy * kappa));
    z = (invAI * ( z * gamma +  vfz * kappa));

    if( collision == 1 || collision == 3 )
	{
		//std::cout<<"tail collision:"<<x(0)<<","<<y(0)<<","<<z(0)<<std::endl;
	   x(0) = tail_index[0];
	   y(0) = tail_index[1];
	   z(0) = tail_index[2];
	}

	if( collision == 2 || collision == 3 )
	{
		//std::cout<<"head collision:"<<x(end)<<","<<y(end)<<","<<z(end)<<std::endl;
	   x(end) = head_index[0];
	   y(end) = head_index[1];
	   z(end) = head_index[2];
	}
	
	for( unsigned int k = 0; k < x.size(); k++ )
	{
	    Cu.Pt[k].x = x(k);
        Cu.Pt[k].y = y(k);
		Cu.Pt[k].z = z(k);
		Cu.Pt[k].check_out_of_range_3D(SM,SN,SZ);
	}

     if( boundary_head == 1 && boundary_tail == 1)
	 {
	   hit_boundary = true;
	   break;
	 }
  

   }
   
   //check for NaN
   for( int i = 0; i < Cu.GetSize(); i++ )
   {
     if( isnan(Cu.Pt[i].x) || isnan(Cu.Pt[i].y) || isnan(Cu.Pt[i].z) )
		 return;
   }

   //resampling
   //Cu.curveinterp_3D((float)pt_distance);
   //resampling
   vnl_vector<float> vnl_Ru(Ru.size());
   for( int i = 0; i < Ru.size(); i++ )
   {
     vnl_Ru(i) = Ru[i];
   }
   vnl_Ru = Cu.curveinterp_4D((float)pt_distance, vnl_Ru);
   Ru.clear();
   for( int i = 0; i < vnl_Ru.size(); i++ )
   {
     Ru.push_back(vnl_Ru(i));
   }


}



void SnakeClass::OpenSnakeDeform(float alpha, int ITER, int pt_distance, float beta, float kappa, float gamma, int N_Active, bool freeze_body)
{

   if( ITER == 0 )
	   return;
   if( collision == 3)
	   return;


   alpha = 0.1;
   beta = 0;
   freeze_body = false;

   typedef itk::VectorLinearInterpolateImageFunction< 
                       GradientImageType, float >  GradientInterpolatorType;
   GradientInterpolatorType::Pointer interpolator = GradientInterpolatorType::New();
   interpolator->SetInputImage(IM->IGVF);

   int SM = IM->I->GetLargestPossibleRegion().GetSize()[0];
   int SN = IM->I->GetLargestPossibleRegion().GetSize()[1];
   int SZ = IM->I->GetLargestPossibleRegion().GetSize()[2];
  
   int N = Cu.GetSize();
   vnl_matrix<float> A = makeOpenA(alpha, beta, N);
   vnl_matrix<float> I(N,N);
   I.set_identity();
   vnl_matrix<float> invAI = vnl_matrix_inverse<float>( A + I * gamma);

   for( int i = 0; i < ITER; i++ )
   {

     vnl_vector<float> vfx(N);
     vnl_vector<float> vfy(N);
	 vnl_vector<float> vfz(N);
     vnl_vector<float> x(N);
     vnl_vector<float> y(N);
	 vnl_vector<float> z(N);
	 vfx.fill(0);
	 vfy.fill(0);
	 vfz.fill(0);

	 for( int j = 0; j < Cu.GetSize(); j++ )
	 {
	    Cu.Pt[j].check_out_of_range_3D(SM,SN,SZ);

		x(j) = Cu.GetPt(j).x;
		y(j) = Cu.GetPt(j).y;
        z(j) = Cu.GetPt(j).z;

        GradientImageType::IndexType index; 
        index[0] = (x(j));
		index[1] = (y(j));
		index[2] = (z(j));
        
       if( freeze_body )
	   {
	   	if( j >= Cu.GetSize() - N_Active || j <= N_Active ) 
	    {
		 GradientPixelType gvf = interpolator->EvaluateAtIndex(index);
         vfx(j) = gvf[0];
         vfy(j) = gvf[1];
	 	 vfz(j) = gvf[2];
	    }
	   }
	   else
	   {
	     GradientPixelType gvf = interpolator->EvaluateAtIndex(index);
         vfx(j) = gvf[0];
         vfy(j) = gvf[1];
	 	 vfz(j) = gvf[2];
	   }
	 }

	 //disable the external froce of ends with collision
	 if( collision == 1 )
	 {
		 vfx(0) = 0;
		 vfy(0) = 0;
		 vfz(0) = 0;
     }
	 else if( collision == 2 )
	 { 
	     vfx(N-1) = 0;
		 vfy(N-1) = 0;
		 vfz(N-1) = 0;
	 }

    x = (invAI * ( x * gamma +  vfx * kappa));
    y = (invAI * ( y * gamma +  vfy * kappa));
    z = (invAI * ( z * gamma +  vfz * kappa));

    for( unsigned int k = 0; k < x.size(); k++ )
	{
	  if( collision == 1 && k == 0 || collision == 2 && k == x.size()-1)
		  continue;

	    Cu.Pt[k].x = x(k);
        Cu.Pt[k].y = y(k);
		Cu.Pt[k].z = z(k);
		Cu.Pt[k].check_out_of_range_3D(SM,SN,SZ);
	}

 
   }
   
   //check for NaN
   for( int i = 0; i < Cu.GetSize(); i++ )
   {
     if( isnan(Cu.Pt[i].x) || isnan(Cu.Pt[i].y) || isnan(Cu.Pt[i].z) )
		 return;
   }

   //resampling
   //Cu.curveinterp_3D((float)pt_distance);

   //resampling
   vnl_vector<float> vnl_Ru(Ru.size());
   for( int i = 0; i < Ru.size(); i++ )
   {
     vnl_Ru(i) = Ru[i];
   }
   vnl_Ru = Cu.curveinterp_4D((float)pt_distance, vnl_Ru);
   Ru.clear();
   for( int i = 0; i < vnl_Ru.size(); i++ )
   {
     Ru.push_back(vnl_Ru(i));
   }

  for( int k = 0; k < Cu.GetSize(); k++ )
  {
		Cu.Pt[k].check_out_of_range_3D(SM,SN,SZ);
  }

}


vnl_matrix<float> SnakeClass::makeOpenA(float alpha, float beta, int N)
{
   vnl_vector<float> Alpha(N);
   Alpha.fill(alpha);
   vnl_vector<float> Beta(N);
   Beta.fill(beta);
   Beta(0) = 0;
   Beta(Beta.size()-1) = 0;

   vnl_matrix<float> A(N,N);
   A.fill(0);
  
   for( int i = 0; i < N; i++ )
   {
     int iplus1 = i + 1;
	 int iplus2 = i + 2;
	 int iminus1 = i - 1;
	 int iminus2 = i - 2;

	 if( iminus1 <= -1 ) {iminus1 = -1 * iminus1;}
	 if( iminus2 <= -1 ) {iminus2 = -1 * iminus2;}
     if( iplus1 > N-1 ) {iplus1 = 2 * (N-1) - iplus1;}
     if( iplus2 > N-1 ) {iplus2 = 2 * (N-1) - iplus2;}

     A(i, iminus2) = A(i, iminus2) + Beta(iminus1);
     A(i, iminus1) = A(i, iminus1) - 2 * Beta(iminus1) - 2 * Beta(i) - Alpha(i);
     A(i, i) = A(i,i) + Alpha(i) + Alpha(iplus1) + 4 * Beta(i) + Beta(iminus1) + Beta(iplus1);
     A(i, iplus1) = A(i, iplus1) - Alpha(iplus1) - 2 * Beta(i) - 2 * Beta(iplus1);
     A(i, iplus2) = A(i, iplus2) + Beta(iplus1);
   }

  return A;
}

void SnakeClass::Grow_Snake_Point()
{
   int SM = IM->V1->GetLargestPossibleRegion().GetSize()[0];
   int SN = IM->V1->GetLargestPossibleRegion().GetSize()[1];
   int SZ = IM->V1->GetLargestPossibleRegion().GetSize()[2];

   int iter_num = 2;

   GradientImageType::IndexType index; 

   if(Cu.NP >= 3)
   {
     for(int i = 1; i < iter_num; i++)
	 {
	   int NP = Cu.NP;
	   //grow head
	   Cu.Pt[NP-1].check_out_of_range_3D(SM,SN,SZ);
	   index[0] = (Cu.Pt[NP-1].x);
       index[1] = (Cu.Pt[NP-1].y);
	   index[2] = (Cu.Pt[NP-1].z);
	   GradientPixelType first_p_direction = IM->V1->GetPixel(index);
	   GradientPixelType flow_direction = IM->IGVF->GetPixel(index);

	   //check direction
	   Vector3D V_V1(first_p_direction[0],first_p_direction[1],first_p_direction[2]);
	   Vector3D V_S;
	   V_S.SetEndPoints(Cu.Pt[NP-2],Cu.Pt[NP-1]);
	   
	   if( V_S.GetDProduct(V_V1) < 0 )
	   {
	      first_p_direction[0] *= -1;
          first_p_direction[1] *= -1;
		  first_p_direction[2] *= -1;
	   }
	   Point3D temp(first_p_direction[0],first_p_direction[1],first_p_direction[2]);
	   Point3D temp0(flow_direction[0],flow_direction[1],flow_direction[2]);
	   temp = Cu.Pt[NP-1] + temp * 2 + temp0;
	   temp.check_out_of_range_3D(SM,SN,SZ);
	   //std::cout<<"head temp:"<<temp.x<<","<<temp.y<<","<<temp.z<<std::endl;
       Cu.AddPt(temp);

       //grow tail
	   Cu.Pt[0].check_out_of_range_3D(SM,SN,SZ);
	   index[0] = (Cu.Pt[0].x);
       index[1] = (Cu.Pt[0].y);
	   index[2] = (Cu.Pt[0].z);
	   first_p_direction = IM->V1->GetPixel(index);
	   flow_direction = IM->IGVF->GetPixel(index);
	   //check direction
	   Vector3D V_V2(first_p_direction[0],first_p_direction[1],first_p_direction[2]);
	   V_S.SetEndPoints(Cu.Pt[1],Cu.Pt[0]);
	   if( V_S.GetDProduct(V_V2) < 0 )
	   {
	      first_p_direction[0] *= -1;
          first_p_direction[1] *= -1;
		  first_p_direction[2] *= -1;
	   }
	   Point3D temp1(first_p_direction[0],first_p_direction[1],first_p_direction[2]);
	   Point3D temp2(flow_direction[0],flow_direction[1],flow_direction[2]);
	   temp1 = Cu.Pt[0] + temp1 + temp2;
	   //std::cout<<"tail temp1:"<<temp1.x<<","<<temp1.y<<","<<temp1.z<<std::endl;
	   temp1.check_out_of_range_3D(SM,SN,SZ);
       Cu.AddTailPt(temp1);
	   //std::cout<<"check point"<<std::endl;
	 }

   }
}


SnakeListClass::SnakeListClass(void)
{
     NSnakes = 0;
	 //only supports 1024 snakes
	 Snakes = new SnakeClass[4000];
	 valid_list.set_size(4000);
	 valid_list.fill(1);
}

void SnakeListClass::SetNSpace(int N)
{
     NSnakes = 0;
	 Snakes = new SnakeClass[N];
	 valid_list.set_size(N);
	 valid_list.fill(1);
}

SnakeListClass SnakeListClass::operator=(SnakeListClass SnakeList)
{ 
	 NSnakes = 0;
	 for( int i = 0; i<SnakeList.NSnakes; i++)
	 {
		 AddSnake(SnakeList.Snakes[i]);
	 }
	 
	 valid_list = SnakeList.valid_list;

	 return *this;
}

void SnakeListClass::RemoveSnake(int idx)
{
	 valid_list(idx) = 0;
     IM->ImCoding( Snakes[idx].Cu, 0, false );
}

void SnakeListClass::SplitSnake(int idx_snake, int idx_pt)
{
     SnakeClass temp;

 	 for( int i = idx_pt + 1; i < Snakes[idx_snake].Cu.GetSize(); i++ )
	 {
	 	temp.Cu.AddPt( Snakes[idx_snake].Cu.Pt[i] );
		temp.Ru.push_back( Snakes[idx_snake].Ru[i] );
	 }

     IM->ImCoding( temp.Cu, NSnakes+1, false );
	 
	 if( temp.Cu.NP >= 2 )
	 {
	  AddSnake(temp);
	 }

	 //Snakes[idx_snake].Cu.NP = idx_pt;
	 Snakes[idx_snake].Cu.Resize(idx_pt);
	 Snakes[idx_snake].Ru.resize(idx_pt);

	 //Snakes[idx_snake].BranchPt.NP = 0;
	 Snakes[idx_snake].BranchPt.RemoveAllPts();
	 
	 if( Snakes[idx_snake].Cu.NP < 2 )
		valid_list(idx_snake) = 0;
	 
}

void SnakeListClass::CreateBranch(int idx1, int idx2)
{
     vnl_vector<float> dist_temp(4);
	 Snakes[idx1].BranchPt.RemoveAllPts();
	 Snakes[idx2].BranchPt.RemoveAllPts();
     
	 vnl_vector<float> dist_temp_h2(Snakes[idx2].Cu.NP);
	 vnl_vector<float> dist_temp_t2(Snakes[idx2].Cu.NP);
	 for( int i = 0; i < Snakes[idx2].Cu.NP; i++ )
	 {
	   dist_temp_h2(i) = Snakes[idx1].Cu.GetLastPt().GetDistTo( Snakes[idx2].Cu.Pt[i] );
	   dist_temp_t2(i) = Snakes[idx1].Cu.GetFirstPt().GetDistTo( Snakes[idx2].Cu.Pt[i] );
	 }
	 dist_temp(0) = dist_temp_h2.min_value();
     dist_temp(1) = dist_temp_t2.min_value();

	 vnl_vector<float> dist_temp_h1(Snakes[idx1].Cu.NP);
	 vnl_vector<float> dist_temp_t1(Snakes[idx1].Cu.NP);
	 for( int i = 0; i < Snakes[idx1].Cu.NP; i++ )
	 {
	   dist_temp_h1(i) = Snakes[idx2].Cu.GetLastPt().GetDistTo( Snakes[idx1].Cu.Pt[i] );
	   dist_temp_t1(i) = Snakes[idx2].Cu.GetFirstPt().GetDistTo( Snakes[idx1].Cu.Pt[i] );
	 }

	 dist_temp(2) = dist_temp_h1.min_value();
     dist_temp(3) = dist_temp_t1.min_value();

	 int idx = dist_temp.arg_min();

	 if( idx == 0 )
	 {
	   Snakes[idx1].Cu.AddPt( Snakes[idx2].Cu.Pt[dist_temp_h2.arg_min()] );
	   Snakes[idx1].Ru.push_back( Snakes[idx2].Ru[dist_temp_h2.arg_min()] );
	   Snakes[idx1].BranchPt.AddPt( Snakes[idx2].Cu.Pt[dist_temp_h2.arg_min()] );
	   Snakes[idx1].Cu.curveinterp_3D(Snakes[idx1].Cu.NP);
	 }
	 else if( idx == 1 )
	 {
	   Snakes[idx1].Cu.AddTailPt( Snakes[idx2].Cu.Pt[dist_temp_t2.arg_min()] );
	   Snakes[idx1].Ru.push_back( Snakes[idx2].Ru[dist_temp_h2.arg_min()] );
	   Snakes[idx1].BranchPt.AddPt( Snakes[idx2].Cu.Pt[dist_temp_t2.arg_min()] );
	   Snakes[idx1].Cu.curveinterp_3D(Snakes[idx1].Cu.NP);
	 }
	 else if( idx == 2 )
	 {
       Snakes[idx2].Cu.AddPt( Snakes[idx1].Cu.Pt[dist_temp_h1.arg_min()] );
	   Snakes[idx2].Ru.push_back( Snakes[idx1].Ru[dist_temp_h1.arg_min()] );
	   Snakes[idx2].BranchPt.AddPt( Snakes[idx1].Cu.Pt[dist_temp_h1.arg_min()] );
	   Snakes[idx2].Cu.curveinterp_3D(Snakes[idx2].Cu.NP);
	 }
	 else
	 {
	   Snakes[idx2].Cu.AddTailPt( Snakes[idx1].Cu.Pt[dist_temp_t1.arg_min()] );
	   Snakes[idx2].Ru.push_back( Snakes[idx1].Ru[dist_temp_h1.arg_min()] );
	   Snakes[idx2].BranchPt.AddPt( Snakes[idx1].Cu.Pt[dist_temp_t1.arg_min()] );
	   Snakes[idx2].Cu.curveinterp_3D(Snakes[idx2].Cu.NP);
	 }
	 
}

void SnakeListClass::MergeSnake(int idx1, int idx2, bool im_coding)
{
     vnl_vector<float> dist_temp(4);

	 //Snakes[idx1].BranchPt.NP = 0;
	 //Snakes[idx2].BranchPt.NP = 0;
	 Snakes[idx1].BranchPt.RemoveAllPts();
	 Snakes[idx2].BranchPt.RemoveAllPts();

	 dist_temp(0) = Snakes[idx1].Cu.GetLastPt().GetDistTo( Snakes[idx2].Cu.GetLastPt() );
	 dist_temp(1) = Snakes[idx1].Cu.GetLastPt().GetDistTo( Snakes[idx2].Cu.GetFirstPt() );
	 dist_temp(2) = Snakes[idx1].Cu.GetFirstPt().GetDistTo( Snakes[idx2].Cu.GetLastPt() );
	 dist_temp(3) = Snakes[idx1].Cu.GetFirstPt().GetDistTo( Snakes[idx2].Cu.GetFirstPt() );

	 int idx = dist_temp.arg_min();

	 

     if( idx == 0 )
	 {
	   //find closest point
	   vnl_vector<float> dist_temp1( Snakes[idx2].Cu.GetSize() );

	   for( int i = 0; i < Snakes[idx2].Cu.GetSize(); i++ )
	   {
		   dist_temp1(i) = Snakes[idx1].Cu.GetLastPt().GetDistTo( Snakes[idx2].Cu.Pt[i] );
	   }
	   int pt_idx = dist_temp1.arg_min();
	   
	   Snakes[idx1].Cu.RemovePt(); //remove the last point

	   //transfer points
	   for( int i = pt_idx; i >= 0; i-- )
	   {
		   Snakes[idx1].Cu.AddPt( Snakes[idx2].Cu.Pt[i] );
		   Snakes[idx1].Ru.push_back( Snakes[idx2].Ru[i] );
	   }
	 }
	 else if( idx == 1 )
	 {
	   //find closest point
	   vnl_vector<float> dist_temp1( Snakes[idx2].Cu.GetSize() );

	   for( int i = 0; i < Snakes[idx2].Cu.GetSize(); i++ )
	   {
		   dist_temp1(i) = Snakes[idx1].Cu.GetLastPt().GetDistTo( Snakes[idx2].Cu.Pt[i] );
	   }
	   int pt_idx = dist_temp1.arg_min();
	   
	   Snakes[idx1].Cu.RemovePt(); //remove the last point

	   //transfer points
       for( int i = pt_idx; i < Snakes[idx2].Cu.GetSize(); i++ )
	   {
		   Snakes[idx1].Cu.AddPt( Snakes[idx2].Cu.Pt[i] );
		   Snakes[idx1].Ru.push_back( Snakes[idx2].Ru[i] );
	   }
	 }
	 else if( idx == 2 )
	 {
	   //find closest point
	   vnl_vector<float> dist_temp1( Snakes[idx2].Cu.GetSize() );

	   for( int i = 0; i < Snakes[idx2].Cu.GetSize(); i++ )
	   {
		   dist_temp1(i) = Snakes[idx1].Cu.GetFirstPt().GetDistTo( Snakes[idx2].Cu.Pt[i] );
	   }
	   int pt_idx = dist_temp1.arg_min();
	   
	   Snakes[idx1].Cu.Flip();
	   Snakes[idx1].Cu.RemovePt(); //remove the last point
       Snakes[idx1].Cu.Flip();

	   //transfer points
       for( int i = pt_idx; i >= 0; i-- )
	   {
		   Snakes[idx1].Cu.AddTailPt( Snakes[idx2].Cu.Pt[i] );
		   //Snakes[idx1].Ru.push_back( Snakes[idx2].Ru[i] );
		   std::vector<float>::iterator it;
		   it = Snakes[idx1].Ru.begin();
           it = Snakes[idx1].Ru.insert ( it , Snakes[idx2].Ru[i] );
	   }
	 }
	 else if( idx == 3 )
	 {
	   //find closest point
	   vnl_vector<float> dist_temp1( Snakes[idx2].Cu.GetSize() );

	   for( int i = 0; i < Snakes[idx2].Cu.GetSize(); i++ )
	   {
		   dist_temp1(i) = Snakes[idx1].Cu.GetFirstPt().GetDistTo( Snakes[idx2].Cu.Pt[i] );
	   }
	   int pt_idx = dist_temp1.arg_min();
	   
	   Snakes[idx1].Cu.Flip();
	   Snakes[idx1].Cu.RemovePt(); //remove the last point
       Snakes[idx1].Cu.Flip();

	   //transfer points
       for( int i = pt_idx; i < Snakes[idx2].Cu.GetSize(); i++ )
	   {
		   Snakes[idx1].Cu.AddTailPt( Snakes[idx2].Cu.Pt[i] );
		   //Snakes[idx1].Ru.push_back( Snakes[idx2].Ru[i] );
		   std::vector<float>::iterator it;
		   it = Snakes[idx1].Ru.begin();
           it = Snakes[idx1].Ru.insert ( it , Snakes[idx2].Ru[i] );
	   }
	 }

	 valid_list(idx2) = 0;

	 if( im_coding )
	 {
	  IM->ImCoding( Snakes[idx1].Cu, idx1+1, false );
	 }

	  //for( int i = 0; i < Snakes[idx1].Cu.GetSize(); i++ )
	  // {
	  //   std::cout<<"Snakes[idx1].Cu:"<<Snakes[idx1].Cu.Pt[i].x<<","<<Snakes[idx1].Cu.Pt[i].y<<","<<Snakes[idx1].Cu.Pt[i].z<<std::endl;
	  //}

}

SnakeClass SnakeListClass::GetSnake(int idx)
{   
     return Snakes[idx];
}

void SnakeListClass::AddSnake(SnakeClass snake)
{
     Snakes[NSnakes] = snake;
	 Snakes[NSnakes].SetImage(IM);
	 NSnakes++;
}

void SnakeListClass::SetImage(ImageOperation *I_Input)
{
     IM = I_Input;
}

bool isnan( float x ) 
{ 
	return (x) != (x); 
};

bool isinf( float x )
{
    return std::numeric_limits<float>::has_infinity &&
           x == std::numeric_limits<float>::infinity() ;
}

float norm_density(float x, float mu, float sigma)
{
   float pi = 3.1415926;
   float p;
   p = 1/MAX(sigma * sqrt(2*pi),std::numeric_limits<float>::epsilon()) * exp( -pow(x-mu,2)/MAX(2*pow(sigma,2),std::numeric_limits<float>::epsilon()));
   return p;
}