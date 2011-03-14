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

#include "OpenSnakeTracer.h"

OpenSnakeTracer::OpenSnakeTracer()
{
   tracing = false;
   tracing_thread = new TracingThread;
   old_snake_length = 0;
}

void OpenSnakeTracer::Init()
{
	SnakeList.NSnakes = 0;
	SnakeList.valid_list.fill(1);
	//SnakeList.SetNSpace(IM->SeedPt.GetSize());
	SnakeList.SetImage(IM);
}

void OpenSnakeTracer::Open_Curve_Snake_Tracing()
{
  if( !options.parallel_tracing )
  {
    //in case of batch tracing, remove the seed convered by snakes growing from seed snakes first
	if( SnakeList.NSnakes != 0 && IM->visit_label.sum() == 0 )
	{
	 for( unsigned int i = 0; i < IM->visit_label.size(); i++ )
     {
      if( IM->visit_label(i) == 0)
	  {
	    int NS = SnakeList.NSnakes;
	    for( int j = 0; j < NS ; j++ )
		{
			if( SnakeList.valid_list(j) == 0 )
				continue;
		   bool removal = false;
           for( int k = 0; k < SnakeList.Snakes[j].Cu.NP ; k++ )
		   {
			if( IM->SeedPt.Pt[i].GetDistTo( SnakeList.Snakes[j].Cu.Pt[k] ) <= options.remove_seed_range )
			{
			  IM->visit_label(i) = 1;
			  removal = true;
			  break;
			}
		   }

		   if( removal )
             break;
		}
	  }
	 }
	}


    //eliminate seed points covered by traced snakes
   if( SnakeList.NSnakes != 0 && old_snake_length != SnakeList.NSnakes )
   {
    for( unsigned int i = 0; i < IM->visit_label.size(); i++ )
    {
      if( IM->visit_label(i) == 0 && SnakeList.NSnakes != 0 )
	  {
	    int NS = SnakeList.NSnakes;
	    for( int j = 0; j < SnakeList.Snakes[ NS - 1 ].Cu.GetSize() ; j++ )
		{
			if( IM->SeedPt.Pt[i].GetDistTo( SnakeList.Snakes[ NS - 1 ].Cu.Pt[j] ) <= options.remove_seed_range )
			{
			  IM->visit_label(i) = 1;
			  break;
			}
		}
	  }
	}

	old_snake_length = SnakeList.NSnakes;

   }
   
   std::cout<<"#Seeds Left:"<< IM->SeedPt.GetSize() - IM->visit_label.sum() << std::endl;

   int visit_label_sum = IM->visit_label.sum();

   current_seeds.RemoveAllPts();

   for( unsigned int i = 0; i < IM->visit_label.size(); i++ )
   {
      if( IM->visit_label(i) == 0 )
	  {
		  current_seeds.AddPt(IM->SeedPt.Pt[i]);
		  IM->current_seed_idx = i;
		  //IM->visit_label(i) = 1;
		  break;
	  }
   }

   if( visit_label_sum != IM->SeedPt.GetSize() )
   {
      Cast_Open_Snake_3D( current_seeds, false );
   }
   else
   {
	   std::cout<<"--------------Processing Finished--------------"<<std::endl;
      tracing_thread->stop();
	  
   }



  }
  else
  {
     
    tracing_thread->setParas( &SnakeList, IM->SeedPt, IM, options);

    tracing = true;

    tracing_thread->start();
	std::cout<<"check point 1"<<std::endl;
  }

}

void OpenSnakeTracer::Cast_Open_Snake_3D( PointList3D seeds, bool manual_seed )
{

   //int iter_num = 50;
  
   //expand 2D seed to 3D
 if( manual_seed )
 {
   //std::cout<<"manual seed"<<std::endl;
   int SM = IM->I->GetLargestPossibleRegion().GetSize()[0];
   int SN = IM->I->GetLargestPossibleRegion().GetSize()[1];
   int SZ = IM->I->GetLargestPossibleRegion().GetSize()[2];

   PointList3D temp_seeds;

  for( int j = 0; j < seeds.NP; j++ )
  {
    
   if( seeds.Pt[j].z != 0 )
   {
	   temp_seeds.AddPt(seeds.Pt[j]);
	   continue;
   }

   ProbImageType::IndexType index; 
   ImageType::IndexType index1;
   index[0] = seeds.Pt[j].x;
   index1[0] = seeds.Pt[j].x;
   index[1] = seeds.Pt[j].y;
   index1[1] = seeds.Pt[j].y;
   vnl_vector<float> saliency(SZ);
   saliency.fill(0);

   bool skeleton_point = false;

     for( int i = 0; i < SZ; i++ )
	 {
       
	   index[2] = i;
	   index1[2] = i;

	   /*if( IM->SBW->GetPixel( index1 ) == 1 )
	   {
		  seeds.Pt[j].z = i;
		  seeds.Pt[j].check_out_of_range_3D(SM,SN,SZ);
		  temp_seeds.AddPt(seeds.Pt[j]);
		  skeleton_point = true;
		  break;
	   }*/
	  
	  //if( IM->VBW->GetPixel( index1 ) == 1 )  
		  saliency(i) = IM->I->GetPixel( index );

	 }
   
     if( !skeleton_point )
	 {
	  int idx = saliency.arg_max();
	  
	  //if( saliency.max_value() == 0 && j > 0 )
	  //{
	  //	seeds.Pt[j].z = seeds.Pt[j-1].z;
      //  seeds.Pt[j].check_out_of_range_3D(SM,SN,SZ);
	  //break;
	  //}
	  seeds.Pt[j].z = idx;
	  seeds.Pt[j].check_out_of_range_3D(SM,SN,SZ);

      //if( saliency.arg_max() != 0 )
		  temp_seeds.AddPt(seeds.Pt[j]);
	 }
  }

    if( temp_seeds.NP == 2 )
	{
		seeds.RemoveAllPts();
		seeds.AddPt(temp_seeds.Pt[0]);
	}
	else
	{
        seeds = temp_seeds;
	}

	 tracing_thread->manual_seed = true;
 }
 else
 {
     tracing_thread->manual_seed = false;
 }
   
   tracing_thread->setParas( &SnakeList, seeds, IM, options);

   //in case of going back to preprocessing step, stopped is ture and hence it needs to be reset
   if( SnakeList.NSnakes == 0 )
   {
	   tracing_thread->stopped = false;
   }

   tracing = true;

   tracing_thread->start();

  //int i = 0;
  //while( i < iter_num )
  //{
 //	 std::cout<<"iteration#"<<i<<std::endl;
 //  snake.OpenSnakeStretch( alpha, ITER, pt_distance, beta, kappa, gamma, stretchingRatio );
 //  i++;
 // }
}

void OpenSnakeTracer::SetImage(ImageOperation *I_Input)
{
   IM = I_Input;
   //visit_label.set_size( IM->SeedPt.GetSize() );
   //visit_label.fill(0);
}

void OpenSnakeTracer::setParas(int pt_distance, float gamma, float stretchingRatio, float minimum_length, 
							   int collision_dist, int remove_seed_range, int deform_iter, bool automatic_merging, 
							   int max_angle, int seed_expansion, bool freeze_body, int s_force, float repeat_ratio, int repeat_dist, int tracing_model, bool parallel_tracing)
{
   options.alpha = 0.1;
   options.stretch_iter = 5;
   options.pt_distance = pt_distance;
   options.beta = 0.1;
   options.kappa = 1;
   options.gamma = gamma;
   options.stretchingRatio = stretchingRatio;
   //options.struggle_dist = 2;
   options.struggle_dist = 0.05;
   options.struggle_th = 15;
   options.minimum_length = minimum_length;
   options.iter_num = 100;
   options.remove_seed_range = remove_seed_range;
   options.deform_iter = deform_iter;
   options.collision_dist = collision_dist;
   options.n_active = 3;
   options.automatic_merging = automatic_merging;
   options.max_angle = max_angle;
   options.seed_expansion = seed_expansion;
   options.freeze_body = freeze_body;
   options.s_force = s_force;
   options.repeat_ratio = repeat_ratio;
   options.repeat_dist = repeat_dist;
   options.tracing_model = tracing_model;
   options.parallel_tracing = parallel_tracing;
}
   
void OpenSnakeTracer::Refine_Branch_Point()
{
	for( int i = 0; i < SnakeList.NSnakes; i++ )
	{
		if( SnakeList.valid_list(i) == 1 )
	      SnakeList.Snakes[i].Branch_Adjustment();
	}

}


TracingThread::TracingThread()
: QThread()
{
   stopped = false;
   ToSuspend = false;
}

void TracingThread::setParas(SnakeListClass *sl, PointList3D s, ImageOperation *im, OptionsStruct io)
{
	Snakes = sl;
	//snake = s;
	seeds = s;
	IM = im;
	options = io;
}

void TracingThread::run()
{

   if( stopped )
	   return;
   if( seeds.NP == 0 )
   {
       emit_traced_signal();
	   return;
   }

//std::cout<<"check point 2"<<std::endl;

if( !options.parallel_tracing )
{
   //SnakeClass *snake;
   //snake = new SnakeClass;
	Snakes->Snakes[Snakes->NSnakes].Set_Seed_Point(seeds);
    Snakes->Snakes[Snakes->NSnakes].BranchPt.RemoveAllPts();
    Snakes->Snakes[Snakes->NSnakes].SetImage(IM);
    Snakes->Snakes[Snakes->NSnakes].SetTracedSnakes(Snakes);
	Snakes->Snakes[Snakes->NSnakes].collision = 0;

   //seed expansion for initializing the tracing
 if( Snakes->Snakes[Snakes->NSnakes].Cu.GetSize() < 3 )
 {
   /*if( options.seed_expansion == 0 )
   {
	  Snakes->Snakes[Snakes->NSnakes].Skeleton_Expansion();
      if( Snakes->Snakes[Snakes->NSnakes].Cu.GetSize() < 3 )
      {
      Snakes->Snakes[Snakes->NSnakes].Set_Seed_Point(seeds);
      Snakes->Snakes[Snakes->NSnakes].Expand_Seed_Point(3);
      }
   } */
   if ( options.seed_expansion == 0 )
   { 
     Snakes->Snakes[Snakes->NSnakes].Expand_Seed_Point(3);
   }
 }
 else
 {
	Snakes->Snakes[Snakes->NSnakes].Cu.curveinterp_3D(options.pt_distance);
	//initialize the radius
    Snakes->Snakes[Snakes->NSnakes].Ru.clear();
    for( int j = 0; j < Snakes->Snakes[Snakes->NSnakes].Cu.GetSize(); j++ )
    {
     Snakes->Snakes[Snakes->NSnakes].Ru.push_back(1);
    }
	//seed snake should be deformed first
    Snakes->Snakes[Snakes->NSnakes].OpenSnakeDeform( options.alpha, 20, options.pt_distance, options.beta, options.kappa, options.gamma, options.n_active, false );
 }

  int i = 0;
  int struggle_label = 0;
  bool invalid = false;

  //Point3D old_head;
  //Point3D old_tail;

  //old_head = Snakes->Snakes[Snakes->NSnakes].Cu.GetLastPt();
  //old_tail = Snakes->Snakes[Snakes->NSnakes].Cu.GetFirstPt();

  //for( int j = 0; j < Snakes->Snakes[Snakes->NSnakes].Cu.GetSize(); j++ )
  //{
  //	  Snakes->Snakes[Snakes->NSnakes].Cu.Pt[j].Print();
  //}

  ImageType::IndexType idx;
  idx[0] = Snakes->Snakes[Snakes->NSnakes].Cu.GetMiddlePt().x;
  idx[1] = Snakes->Snakes[Snakes->NSnakes].Cu.GetMiddlePt().y;
  idx[2] = Snakes->Snakes[Snakes->NSnakes].Cu.GetMiddlePt().z;
  int snake_id = IM->VBW->GetPixel(idx);


  //initialize the radius
  Snakes->Snakes[Snakes->NSnakes].Ru.clear();
  for( int j = 0; j < Snakes->Snakes[Snakes->NSnakes].Cu.GetSize(); j++ )
  {
    Snakes->Snakes[Snakes->NSnakes].Ru.push_back(1);
  }
   
  if( options.tracing_model == 2 || options.tracing_model == 3 )
  {
     //Snakes->Snakes[Snakes->NSnakes].OpenSnake_Init_4D(options.alpha, options.stretch_iter, options.beta, options.kappa, options.gamma, options.pt_distance);
     //IM->ImComputeBackgroundModel();
	 //IM->ImComputeForegroundModel(Snakes->Snakes[Snakes->NSnakes].Cu);
  }

  while( i < options.iter_num && struggle_label <= options.struggle_th && !Snakes->Snakes[Snakes->NSnakes].hit_boundary )
  {
	  float old_dist = Snakes->Snakes[Snakes->NSnakes].Cu.GetLength();

	  Snakes->Snakes[Snakes->NSnakes].OpenSnakeDeform( options.alpha, options.deform_iter, options.pt_distance, options.beta, options.kappa, options.gamma, options.n_active, options.freeze_body );

	  if( options.tracing_model == 0 )
	  {
	   Snakes->Snakes[Snakes->NSnakes].OpenSnakeStretch( options.alpha, options.stretch_iter, options.pt_distance, options.beta, options.kappa, options.gamma,
	  	  options.stretchingRatio, options.collision_dist, options.n_active, options.minimum_length, options.automatic_merging, options.max_angle, options.freeze_body, options.s_force, 0);
	  }
	  else
	  {
	   Snakes->Snakes[Snakes->NSnakes].OpenSnakeStretch_4D( options.alpha, options.stretch_iter, options.pt_distance, options.beta, options.kappa, options.gamma,
	  	  options.stretchingRatio, options.collision_dist, options.minimum_length, options.automatic_merging, options.max_angle, options.freeze_body, options.s_force, 0, options.tracing_model);
	  }

	  float new_dist = Snakes->Snakes[Snakes->NSnakes].Cu.GetLength();

	  if( new_dist > old_dist * ( 1 - options.struggle_dist ) && new_dist < old_dist * ( 1 + options.struggle_dist ) )
	     struggle_label++;
	  else
	     struggle_label = 0;

	  //if( Snakes->Snakes[Snakes->NSnakes].Cu.GetFirstPt().GetDistTo(old_tail) <= options.struggle_dist &&
	  //	  Snakes->Snakes[Snakes->NSnakes].Cu.GetLastPt().GetDistTo(old_head) <= options.struggle_dist)
	  //	  struggle_label++;

	  //if( i > options.struggle_th && (Snakes->Snakes[Snakes->NSnakes].Cu.GetSize() == 3 || Snakes->Snakes[Snakes->NSnakes].Cu.GetLength() <= options.minimum_length) ) 
	  //{
	  //   std::cout<<"invalid seed"<<std::endl;
	  //   invalid = true;
	  //   break;
	  //}

      
	  if( Snakes->Snakes[Snakes->NSnakes].Cu.GetLength() > options.minimum_length )
	  {
		  IM->ImCoding( Snakes->Snakes[Snakes->NSnakes].Cu, 1, true );
	  }

	  i++;
	  
	   emit stretched(Snakes->Snakes[Snakes->NSnakes]);

      //wait for displaying image to complete
      //suspend();

      if(ToSuspend) 
	  {
	 	mutex.lock();
		condition.wait(&mutex);
		mutex.unlock();
      }
	  //msleep(5);
	  //resume();
  }

  if( !invalid )
  {
    if( manual_seed )
	//if( 0 )
	{
	   if( Snakes->Snakes[Snakes->NSnakes].Cu.GetLength() >= options.minimum_length )
	   {
	    //Snakes->AddSnake(*snake);
	    IM->ImCoding( Snakes->Snakes[Snakes->NSnakes].Cu, Snakes->NSnakes+1, false );
		//adjust the branch point
		//Snakes->Snakes[Snakes->NSnakes].Branch_Adjustment();
		Snakes->Snakes[Snakes->NSnakes].Nail_Branch();

        Snakes->NSnakes++;
		IM->ImRefresh_TracingImage();
	   }
	   else
	   {
		   std::cout<<"length than minimum length"<<std::endl;
	   }
	}
	else
	{
	  if( Snakes->Snakes[Snakes->NSnakes].Check_Validity( options.minimum_length, options.repeat_ratio, options.repeat_dist, 0) )
	  {
	   //Snakes->AddSnake(*snake);
	   IM->ImCoding( Snakes->Snakes[Snakes->NSnakes].Cu, Snakes->NSnakes+1, false );
	   //adjust the branch point
	   //Snakes->Snakes[Snakes->NSnakes].Branch_Adjustment();
	   Snakes->Snakes[Snakes->NSnakes].Nail_Branch();
       Snakes->NSnakes++;
	   IM->ImRefresh_TracingImage();
	  }
	}
  }

   if( !manual_seed)
    IM->visit_label[IM->current_seed_idx] = 1;

    //emit snakeTraced(Snakes);
	emit_traced_signal();
}
else //Parallel Tracing
{
    Snakes->SetNSpace(seeds.NP);
	Snakes->NSnakes = seeds.NP;

  for( int i = 0; i < seeds.NP; i++ )
  {
     Snakes->Snakes[i].Set_Seed_Point(seeds.Pt[i]);
	 Snakes->Snakes[i].BranchPt.RemoveAllPts();
     Snakes->Snakes[i].SetImage(IM);
     Snakes->Snakes[i].SetTracedSnakes(Snakes);
	 Snakes->Snakes[i].collision = 0;

     Snakes->Snakes[i].Expand_Seed_Point(3);

	 Snakes->Snakes[i].Ru.clear();
     for( int j = 0; j < Snakes->Snakes[i].Cu.GetSize(); j++ )
    {
      Snakes->Snakes[i].Ru.push_back(1);
    }
  }

//std::cout<<"check point 3"<<std::endl;

  int k = 0;
  while( k < 50 )
  {
    for( int i = 0; i < Snakes->NSnakes; i++ )
	{
	 if( Snakes->valid_list(i) == 0 )
		 continue;

     Snakes->Snakes[i].OpenSnakeDeform( options.alpha, 1, options.pt_distance, options.beta, options.kappa, options.gamma, options.n_active, options.freeze_body );
     Snakes->Snakes[i].OpenSnakeStretch( options.alpha, 1, options.pt_distance, options.beta, options.kappa, options.gamma,
	  	  options.stretchingRatio, options.collision_dist, options.n_active, options.minimum_length, options.automatic_merging, options.max_angle, options.freeze_body, options.s_force, i+1);
	  
	  if( Snakes->Snakes[i].Cu.GetLength() > options.minimum_length )
	  {
		  IM->ImCoding( Snakes->Snakes[i].Cu, i+1 , false );
	  }
	  else
	  {
	    Snakes->valid_list(i) = 0;
	  }
	}

	k++;

	 if(ToSuspend) 
	 {
	 	mutex.lock();
		condition.wait(&mutex);
		mutex.unlock();
     }
	 emit_traced_signal();
  }

  for( int i = 0; i < Snakes->NSnakes; i++ )
  {
     //if( !Snakes->Snakes[i].Check_Validity( options.minimum_length, options.repeat_ratio, options.repeat_dist, i+1 ) )
	 //	 Snakes->valid_list(i) = 0;
	  std::cout<<"snakes_valid_list:"<<Snakes->valid_list(i)<<std::endl;
	  
  }
  emit_traced_signal();
  std::cout<<"snakes:"<<Snakes->NSnakes<<std::endl;
  
}

}

void TracingThread::emit_traced_signal()
{
   if( manual_seed )
   {
    //count the length difference, and use it as the false negative
	 float length1 = 0.0;
	 float length2 = 0.0;
	for(int i = 0; i < Snakes->NSnakes; i++)
	{
	    if( Snakes->valid_list(i) == 0 )
			continue;

		float length = Snakes->Snakes[i].Cu.GetLength();

	   if( i != Snakes->NSnakes-1 )
		length1 += length;

	   length2 += length;
	}

	float fn = length2 - length1;

	emit snakeTraced_manual_seed(fn);
   }
   else
	emit snakeTraced();
}

void TracingThread::suspend()
{
   ToSuspend = true;
}

void TracingThread::resume()
{
  if (!ToSuspend)
     return;
   ToSuspend = false;
   condition.wakeOne();
}

void TracingThread::stop()
{
   if( ToSuspend )
	  resume();
   stopped = true;
}
