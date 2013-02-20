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

void OpenSnakeTracer::Init()
{
	SnakeList.SetNSpace(IM->SeedPt.GetSize());
	SnakeList.SetImage(IM);
	seed_num = 0;
}

void OpenSnakeTracer::Open_Curve_Snake_Tracing()
{

   for( int i = 0; i < IM->SeedPt.GetSize(); i++)
   {
	   std::cout<<"#Seeds Left:"<< IM->SeedPt.GetSize() - i << std::endl;
	   Cast_Open_Snake_3D();
   }
   
}


void OpenSnakeTracer::Cast_Open_Snake_3D()
{

   float alpha = 0; 
   int iter_num = 50;
   int ITER = 5;

   int pt_distance = 2;

   float beta = 0.05;
   float kappa = 1;
   float gamma = 1;
   float stretchingRatio = 3;

   Point3D seed = IM->SeedPt.Pt[seed_num];
   seed_num++;
   std::cout<<"#Seeds Left:"<< IM->SeedPt.GetSize() - seed_num << std::endl;

   SnakeClass snake;
   snake.Set_Seed_Point(seed);
   snake.SetImage(IM);
   snake.Expand_Seed_Point(3);

  int i = 0;
  while( i < iter_num )
  {
	  std::cout<<"iteration#"<<i<<std::endl;
	  snake.OpenSnakeStretch( alpha, ITER, pt_distance, beta, kappa, gamma, stretchingRatio );
	  i++;
  }
   
   SnakeList.AddSnake(snake);


}

void OpenSnakeTracer::SetImage(ImageOperation *I_Input)
{
   IM = I_Input;
}
   	


TracingThread::TracingThread()
: QThread()
{

}

void TracingThread::setSnake(SnakeClass *s)
{
	snake = s;
}

void TracingThread::setIterNum(int num)
{
	iter_num = num;
}

void TracingThread::run()
{

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
   stopped = true;
}
