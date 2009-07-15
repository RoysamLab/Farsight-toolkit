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

//
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <list>
#include <queue>

#include "CONSTANTS.h"
#include "Mytypes.h"
#include "Config.h"
#include "Cpoints.h"
#include "Dllist.h"
#include "Cimage.h"
#include "C3dimage.h"
#include "Cvector.h"
#include "Cvessel.h"
#include "Template.h"
#include "Extern.h"
#include "Soma.h"

extern int giParam_Trace_MaxNumOfConsecutiveZeroResponses;
extern int giParam_Trace_MaxStepSize;

double compute_Q(C3DImage * VesselnessImage, double alpha)
{
	int number_of_traced_points = 0;
	int number_of_segments = 0;
	
	for (int v = 0; v < gTheVessels.m_iNumOfElements; v++)
	{
		if ( gTheVessels.m_apData[v] )
		{
			CLNode<CPoint>* temp = gTheVessels.m_apData[v]->m_Center.head;
			while (temp)
			{
				number_of_traced_points++;
				temp = temp->after;
			}
			number_of_segments++;
		}
	}
	int number_of_chained_pixels = number_of_traced_points - number_of_segments;

	std::vector<double> F_DL_LUT(256);
	std::vector<double> B_DL_LUT(256);
	std::string tmp;

	std::ifstream read_F_DL_LUT("F_DL_LUT.txt");
	if(!read_F_DL_LUT) // use the default hypothesized (inverted \__ __/) exponential distributions 
	{
		double fg_rate = 1.0 / 128.0;
		for (int i = 1; i <= 255; i++)
			F_DL_LUT[i] =  ceil(- ( log((exp(-fg_rate * (255 - i - 1.0) ) - exp( -fg_rate * (255 - i) ))) / log(2.0))); 
		F_DL_LUT[0] = F_DL_LUT[1];
	}
	else // read in the description lengths for the foreground region
	{
		for (int i = 0; i <= 255; i++)
			read_F_DL_LUT >> F_DL_LUT[i];
	}

	std::ifstream read_B_DL_LUT("B_DL_LUT.txt");
	if(!read_B_DL_LUT) // use the default hypothesized exponential distribution
	{
		double bg_rate = 1.0 / 0.78534;
		for (int i = 1; i <= 255; i++)
			B_DL_LUT[i] = ceil(- ( log((exp(-bg_rate * (i - 1.0) ) - exp( -bg_rate * i ))) / log(2.0))); 
		B_DL_LUT[0] = B_DL_LUT[1];
	}
	else // read in the description lengths for the background region
	{
		for (int i = 0; i <= 255; i++)
			read_B_DL_LUT >> B_DL_LUT[i];
	}

	// now get the coverage term
	double term1 = 0.0;
	double p_x_given_m = 0.0;
	//double const_term = 1 / sqrt(2*3.142*1.0);
	double sample_value = 0.0;
	//double asymptotic_dl = 0.0;
	//double min_prob = 1.0 - (254.0/255.0);
	//double max_prob = 254.0 / 255.0;

	int iPadding = gConfig.GetImagePadding();

	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;

	for (int s = iPadding; s < iSlices - iPadding; s++)
	{
		for (int r = 0; r < iRows; r++)
		{
			for (int c = 0; c < iCols; c++)
			{
				sample_value = VesselnessImage->data[s][r][c];
				if( TracedImage[s][r][c] == 0 )
				{	
					p_x_given_m = B_DL_LUT[static_cast<int>(sample_value)];
				}
				else
				{
					p_x_given_m = F_DL_LUT[static_cast<int>(sample_value)];
				}	
				term1 += p_x_given_m;
			}
		}
	}
	
	std::cout << "term1: " << term1 << std::endl;

	// calculate the model encoding description length;

	// 3. calculate the description length
	double cost_per_trace_pixel = ceil( log(27.0) / log(2.0) );
	double cost_per_start_point = 96.0;
	double cost_of_width = 64.0;

	double model_DL = cost_per_trace_pixel * static_cast<double>(number_of_chained_pixels) 
		+ cost_per_start_point * static_cast<double>(number_of_segments)
		+ cost_of_width * static_cast<double>(number_of_traced_points);

	std::cout << "model_DL: " << model_DL << std::endl;
	std::cout << "total DL: " << model_DL + term1 << std::endl;

	double Q;
	if ( Round(model_DL) != 0 )
		Q = (alpha * term1) + ((1.0-alpha) * model_DL);
	else
		Q = (alpha * term1) + ((1.0-alpha) * term1);
	// temp - open result file, read in and write out the column headers if it doesn't exist

	int RSD = gConfig.GetRelativeShiftDistance();	
	int FL = gConfig.GetMinimumTemplateLength();
	int TL = gConfig.GetMaximumTemplateLength();	
	int DF = gConfig.GetDirectionalDegreeOfFreedom();
	int STEP = gConfig.GetMaximumStepSize();
	int CTM = gConfig.GetContrastThresholdMultiplier();
	int MASV = gConfig.GetMaximumAllowedStoppingViolations();
	int GS = gConfig.GetGridSpacing();

	int DL = Round(Q);

	std::string output_path = gConfig.GetOutputPath();
	std::string image_name = gConfig.GetImageName();
	std::string result_filename = output_path + image_name + "_DL_result_appended.txt";
	std::ifstream in(result_filename.c_str());
	if( !in ) // write out headers and current result
	{
		std::ofstream out(result_filename.c_str());
		out << "RSD\tFL\tTL\tDF\tSTEP\tCTM\tMASV\tGS\tDL1\tDL2\tDL\n";
		out << RSD << "\t" 
			<< FL << "\t" 
			<< TL << "\t" 
			<< DF << "\t"
			<< STEP << "\t" 
			<< CTM << "\t" 
			<< MASV << "\t" 
			<< GS << "\t" 
			<< term1 << "\t"
			<< model_DL << "\t"
			<< DL << "\n";
	}
	else	// append
	{
		std::ofstream out(result_filename.c_str(), std::ofstream::app);
		out << RSD << "\t" 
			<< FL << "\t" 
			<< TL << "\t" 
			<< DF << "\t"
			<< STEP << "\t" 
			<< CTM << "\t" 
			<< MASV << "\t" 
			<< GS << "\t" 
			<< term1 << "\t"
			<< model_DL << "\t"
			<< DL << "\n";
	}	

	return Q;
}
