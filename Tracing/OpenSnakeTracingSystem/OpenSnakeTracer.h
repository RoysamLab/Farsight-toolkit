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
  Program:   Open Snake Tracing System
  Autohr:    Yu Wang
  Email: wangy15@rpi.edu
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $
=========================================================================*/

#ifndef OPENSNAKETRACER_H
#define OPENSNAKETRACER_H

#include "itkNumericTraits.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <QTextEdit>
#include "TracingCore/SnakeOperation.h"
#include "qthread.h"
#include "qmutex.h"
#include "qwaitcondition.h"

struct OptionsStruct
{
   float alpha; 
   int stretch_iter;
   int pt_distance;
   float beta;
   float kappa;
   float gamma;
   float stretchingRatio;
   float struggle_dist;
   float struggle_th;
   float minimum_length;
   float iter_num;
   float remove_seed_range;
   float deform_iter;
   int collision_dist;
   int n_active;
   bool automatic_merging;
   int max_angle;
   int seed_expansion;
   bool jumping_gaps;
   int jumping_dist;
   int jumping_angle;
   bool jumping_crossover;
   int crossover_dist;
   bool freeze_body;
   int s_force;
   float repeat_ratio;
   int repeat_dist;
   int tracing_model;
   bool parallel_tracing;
};

class TracingThread : public QThread
{
 Q_OBJECT
 public:

	TracingThread();
    
	void setParas(SnakeListClass *sl, PointList3D s, ImageOperation *im, OptionsStruct io);

    void run();
	
    volatile bool stopped;
	bool ToSuspend;

	PointList3D seeds;
	SnakeListClass *Snakes;
	
	ImageOperation *IM;
    OptionsStruct options;

	QWaitCondition condition;
	QMutex mutex;

	bool manual_seed;
    bool parallel_tracing;

 public slots:
	void suspend();
	void stop();
    void resume();
	void emit_traced_signal();

 signals:
    void stretched(SnakeClass s);
	void snakeTraced(SnakeListClass *s);
	void snakeTraced();
	void snakeTraced_manual_seed(float fn);

 private:

};

class OpenSnakeTracer
{
public:
   
	OpenSnakeTracer();
    SnakeListClass SnakeList;
	
	void Init();
	void setParas(int pt_distance, float gamma, float stretchingRatio, float minimum_length, int collision_dist, 
		          int remove_seed_range, int deform_iter, bool automatic_merging, int max_angle, int seed_expansion,
				  bool freeze_body, int s_force, float repeat_ratio, int repeat_dist, int tracing_model, bool parallel_tracing);
	void SetImage(ImageOperation *I_Input);
	void Open_Curve_Snake_Tracing();
	void Cast_Open_Snake_3D(PointList3D seeds, bool manual_seed);

	void Refine_Branch_Point();
	//Functions for GUI
	QTextEdit *log_viewer;

    OptionsStruct options;

	int old_snake_length;
	bool tracing;

	PointList3D current_seeds;

	ImageOperation *IM;
	//thread for controlling snake stretching and dynamic displaying
	TracingThread *tracing_thread;
};


#endif
