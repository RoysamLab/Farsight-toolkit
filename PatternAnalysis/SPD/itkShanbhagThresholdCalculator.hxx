/*=========================================================================
*
*  Copyright Insight Software Consortium
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*         http://www.apache.org/licenses/LICENSE-2.0.txt
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*=========================================================================*/

#ifndef __itkShanbhagThresholdCalculator_hxx
#define __itkShanbhagThresholdCalculator_hxx

#include "itkShanbhagThresholdCalculator.h"
#include "itkProgressReporter.h"
#include "vnl/vnl_math.h"
#include <fstream>

namespace itk
{
	/*
	* Compute the Shanbhag's threshold
	*/
	template<class THistogram, class TOutput>
	void
		ShanbhagThresholdCalculator<THistogram, TOutput>
		::GenerateData(void)
	{
		const HistogramType * histogram = this->GetInput();
		// histogram->Print(std::cout);
		if ( histogram->GetTotalFrequency() == 0 )
		{
			itkExceptionMacro(<< "Histogram is empty");
		}
		ProgressReporter progress(this, 0, histogram->GetSize(0) );
		if( histogram->GetSize(0) == 1 )
		{
			this->GetOutput()->Set( static_cast<OutputType>( histogram->GetMeasurement(0,0) ) );
		}

		unsigned int size = histogram->GetSize(0);

		const double tolerance = 2.220446049250313E-16;
		int threshold;
		int ih, it;
		int first_bin;
		int last_bin;
		double term;
		double tot_ent;  /* total entropy */
		double min_ent;  /* max entropy */
		double ent_back; /* entropy of the background pixels at a given threshold */
		double ent_obj;  /* entropy of the object pixels at a given threshold */
		std::vector<double> norm_histo(size); /* normalized histogram */
		std::vector<double> P1(size); /* cumulative normalized histogram */
		std::vector<double> P2(size);

		std::ofstream ofs("Shanbhag.txt");
		std::vector<double> smoothedHist(size);
		int total = 0;
		for( int i = 0; i<size; i++)
		{ 
			smoothedHist[i] = static_cast< double >( histogram->GetFrequency(i, 0) );
			total +=  smoothedHist[i];
			ofs <<  smoothedHist[i] <<"\t";
		}
		ofs<< std::endl;

		total = total - smoothedHist[0] ;
		smoothedHist[0] = 0;


		for (ih = 0; (unsigned)ih < size; ih++ )
		{
			norm_histo[ih] = (double)smoothedHist[ih]/total;
			ofs<<  norm_histo[ih]<<"\t";
		}
		ofs<< std::endl;

		for (int i = 0; i<size; i++)
		{
			ofs<< static_cast<double>(histogram->GetMeasurement(i, 0 ) )/ 255<<"\t";
		}
		ofs<< std::endl;


		P1[0]=norm_histo[0];
		P2[0]=1.0-P1[0];
		for (ih = 1; (unsigned)ih < size; ih++ )
		{
			P1[ih]= P1[ih-1] + norm_histo[ih];
			P2[ih]= 1.0 - P1[ih];
		}

		/* Determine the first non-zero bin */
		first_bin=0;
		for (ih = 0; (unsigned)ih < size; ih++ )
		{
			if ( !(vcl_abs(P1[ih])<tolerance))
			{
				first_bin = ih;
				break;
			}
		}

		/* Determine the last non-zero bin */
		last_bin=size - 1;
		for (ih = size - 1; ih >= first_bin; ih-- )
		{
			if ( !(vcl_abs(P2[ih])<tolerance))
			{
				last_bin = ih;
				break;
			}
		}
		std::cout<< first_bin<<"\t"<< last_bin<<std::endl;

		// Calculate the total entropy each gray-level
		// and find the threshold that maximizes it
		threshold =-1;
		min_ent = 0.01; //itk::NumericTraits<double>::max();

		for ( it = last_bin; it >= first_bin; it-- )
		{
			/* Entropy of the background pixels */
			ent_back = 0.0;
			term = 0.5 / P1[it];
			for ( ih = 1; ih <= it; ih++ )
			{ //0+1?
				ent_back -= norm_histo[ih] / P1[it] * vcl_log ( 0.5 + term * (P1[it] - P1[ih - 1] ));
			}

			/* Entropy of the object pixels */
			ent_obj = 0.0;
			term = 0.5 / P2[it];
			for ( ih = it + 1; (unsigned)ih < size; ih++ )
			{
				ent_obj -= norm_histo[ih] / P2[it] * vcl_log ( 0.5 + term * (P2[it] - P2[ih] ));
			}

			/* Total entropy */
			tot_ent = vcl_abs ( ent_back - ent_obj );
			ofs<< tot_ent << "\t";

			if ( tot_ent < min_ent )
			{
				min_ent = tot_ent;
				threshold = it;
				break;
			}
		}

		ofs.close();

		std::cout<< "Threshold bin: "<< threshold<<std::endl;
		this->GetOutput()->Set( static_cast<OutputType>( histogram->GetMeasurement( threshold, 0 ) ) );
	}

} // end namespace itk

#endif
