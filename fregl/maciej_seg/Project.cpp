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

/** @file Project.cpp
*   @brief main program for cell segmentation
*
*   @author Maciej Wotjon
*/
//#include <crtdbg.h>

#include <vector>

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>	//time execution
#include <vnl/algo/vnl_determinant.h>
#include <stdlib.h>	//for wriring/reading data
#include <stdlib.h>	//printf
#include <sstream>	//stringstream

#include "filter.h"
#include "cell_mgr.h"

#include <string>

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"

#include "itkImageAdaptor.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"

//#include "Stackwalker.h"	//mem leak

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

/********************************************************************************/
/*Desc- main program for running nuclear seg, the program takes 19 inputs args. The inputs
/*are described below. The inputs are sets of consecitive numered images that are the slices
/*of cells. The images are loaded, filtered, using ITK, segmented, then merged and a final score is 
/*computd for the segmentation. There are options to seg the image form data already available or 
/*auto train.
/*
/*	19 input arguments
/*	1 - path - path to folder containing images
/*	2 - filename - name of file name of the image stacks, only the name without the numbering
/*	3 - type - type of image, ei extension tif
/*	4 - sstack - starting number of images ie image_nameXX.tif XX->sstack
/*	5 - estack - ending stack number
/*	6 - swidth - width of number ie image_nameXX.tif swidth -> 2 , image_nameXXX.tif swidth -> 3
/*	7 - pT - threshold parameter, gets divided by 1000, int needed for RRS
/*	8 - g - grid size parameter
/*	9 - rM - radius of median filter
/*	10 - rG - radius of morphological operators
/*	11 - min - should always be 0	(param for watershed) (obselete)
/*	12 - max - should always be 254 (param for watershed) (obselete)
/*	13 - step - param for step size of watershed, should be 1 (obselete)
/*	14 - gliadat - path and filename for glia train data
/*	15 - neurdat - path and filename for neuron train data
/*	16 - channel - channel of nuclei stain R-0 G-1 B-2
/*	17 - mode - 0 -> train from data, 1 -> autotrain
/*	18 - loadimages - 1 -> load preprocessed images, 0-> filter images, always filter images, preprocessed sometimes are incorrect
/*	19 - saveoutput - output segmented image->1 , no output->0
/********************************************************************************/
using namespace itk;

int main( int argc , char * argv[] )
{
	if( argc != 17 )
	{
		std::cerr << argc << " Wrong command line arguments" << std::endl;
		return 0;
	}
	int g = 1 ,		//amount of grids
		rM = 1 ,	//median radius
		rG = 1 ,	//morph radius
		min = 0 , 
		max = 254 , 
		step = 1;
	double  pT = 161;
	std::string temp = "_filt.tif";
	
	int sstack = 0, //starting stack
	    estack = 1,	//end stack
	    swidth = 2,	//width of numbers in filename
		channel = 0,	//channel of nuclei stain R-0 G-1 B-2
		mode = 0,	//0 -> train from data, 1 -> autotrain
		loadimages = 0,	//1 -> load preprocessed images, 0-> filter images
		saveoutput = 0;	//output segmented image->1 , no output->0
	
	std::string inputs , path , filename , type, gliadat, neurdat;

	std::ostringstream oss;
	
	path = argv[1];
	
	filename = argv[2];
	
	type = argv[3];
	
	inputs = argv[4];
	std::stringstream( inputs ) >> sstack;
	
	inputs = argv[5];
	std::stringstream( inputs ) >> estack;
	
	inputs = argv[6];
	std::stringstream( inputs ) >> swidth;
	
	inputs = argv[7];
	std::stringstream( inputs ) >> pT;
	
	inputs = argv[8];
	std::stringstream( inputs ) >> g;
	
	inputs = argv[9];
	std::stringstream( inputs ) >> rM;
	
	inputs = argv[10];
	std::stringstream( inputs ) >> rG;
	
        /*
	inputs = argv[11];
	std::stringstream( inputs ) >> min;
	
	inputs = argv[12];
	std::stringstream( inputs ) >> max;
	
	inputs = argv[13];
	std::stringstream( inputs ) >> step;
        */

	gliadat = argv[11];

	neurdat = argv[12];

	inputs = argv[13];
	std::stringstream( inputs ) >> channel;

	inputs = argv[14];
	std::stringstream( inputs ) >> mode;	//0 -> train from data, 1 -> autotrain

	inputs = argv[15];
	std::stringstream( inputs ) >> loadimages;

	inputs = argv[16];
	std::stringstream( inputs ) >> saveoutput;

	if( !(mode == 0 || mode == 1) )
	{
		std::cerr << mode << " Wrong mode" << std::endl;
		return 0;
	}

	//norm threshold param
	pT /= 1000.0f;

	oss << '_' << std::setw(2) << std::setfill('0') << g << '_' << std::setw(2) << std::setfill('0') << rM << '_' << std::setw(2) << std::setfill('0') << rG << '_' << std::setw(6) << std::setfill('0') << pT;
	
	std::cout << "Starting filter " << g << ' ' << rM << ' ' << pT << std::endl;
	filter f(g,g,rM,rG,pT,0.5,channel);	//setup filter to prefilter image

	if( loadimages == 1 )
	{
		std::cout << "Loading Images" << std::endl;
		f.load_images(std::string( "_filt" + oss.str() +".tif"),std::string( "_gradient" + oss.str() +".tif"),std::string( "_intensity" + oss.str() +".tif"),path,filename,"tif",sstack,estack,swidth,loadimages);
	}
	if( loadimages == 0 )
	{
		std::cout << "Filtering Images" << std::endl;
		f.run_filter(filename,path,type,sstack,estack,swidth);		//run filter with path to images and size of stack
		//never do this
		if( saveoutput == 10 )
		{

			temp = "_filt" + oss.str() + ".tif";
			f.save_image(f.getfilt(),temp);

			temp = "_intensity" + oss.str() + ".tif";
			f.save_image(f.getintensity(),temp);

			temp = "_gradient" + oss.str() + ".tif";
			f.save_image(f.getgradient(),temp);
		}
	}

	std::cout << "Starting cell" << std::endl;
	//start cell
	cell_mgr c(f.getsize(),f.getfilt(),f.getintensity(),f.getgradient());

	c.setup_cells(min,max,step,oss.str());	//create cells from labeled image

	if( mode == 0 )
	{
		std::cout << "Getting train data" << std::endl;
		c.get_data(gliadat,neurdat);
	}
	else
	{
		std::cout << "Calc train data" << std::endl;
		c.calc_train_data( 0.95 , 1.0 , 0.9 , 0.2 , 1000 , 8000);

		//merge traning cells
		std::cout << "Build train rag" << std::endl;
		c.build_rag_train();	

		//get rid of cells touching border
		c.remove_border_cells();

		std::cout << "Classifying" << std::endl;
		c.classify_cells();	
		
		std::cout << "Calculating train data" << std::endl;
		//save data into matrix
		c.calc_ng_data();

		//only keep good cells as training cells
		c.calc_mean_vect(2.0/3.0);

		//save good cells into matrix
		c.calc_ng_data();

		//reset cell class labels
		std::cout << "Relabeling Orig Cells" << std::endl;
		c.relabel_orig_cells();		

	}

	//merge and classify cells from train data
	std::cout << "Building 2nd rag" << std::endl;
	c.build_rag(mode,oss.str());	//merge cells

	std::cout << "Removing border Cells" << std::endl;
	c.remove_border_cells();

	if( mode == 1)
	{
		std::cout << "Cleaning up data" << std::endl;
		//c.save_train_data();	
	}	

	c.calc_ng_data();

	//get rid of cells not touching median plane
	c.remove_non_median_cells();

	//correct any mis classfied cells
	std::cout << "Correct Errors" << std::endl;
	c.build_rag(mode,oss.str(),2);

	//relabel cells with consecutive label numbers
	std::cout << "Relabeling" << std::endl;
	c.relabel_cells(true);
	c.calc_ng_data();

	//score the segmentation
	std::cout << "Scoring" << std::endl;
	c.score_foreground_mah( oss.str() );

	c.check_all_cells(0,100,100,1,200);
	if( saveoutput == 1 )
	{
		std::cout << "Labeling and Saving" << std::endl;
		//label border and number cells
		f.label_bound(c.get_cells());
		f.label_image(c.get_cells());
		//save image
		//f.save_image_rgb(std::string("_out_" + oss.str() + ".tif"));
		f.save_image_rgb(std::string(filename + "_out_.tif"));
		//save cell features
		c.output_cells(std::string( filename + ".txt" ));
		//project cells
		std::cout << "Projecting cells" << std::endl;
		c.project_cells();
		//save projected image
		std::cout << "Saving projected" << std::endl;
		f.project_image( filename + "_proj_.tif",c.get_cells());
	}
	return 0;
}
