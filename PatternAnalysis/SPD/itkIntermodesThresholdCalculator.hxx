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
#ifndef __itkIntermodesThresholdCalculator_hxx
#define __itkIntermodesThresholdCalculator_hxx

#include "itkIntermodesThresholdCalculator.h"
#include "itkProgressReporter.h"
#include <fstream>

#define PEAK_NUM 2
namespace itk
{

	template<class THistogram, class TOutput>
	bool
		IntermodesThresholdCalculator<THistogram, TOutput>
		::BimodalTest(const std::vector<double> & h)
	{
		int modes = 0;
		const size_t len = h.size();
		for (size_t k = 1; k < len - 1; k++)
		{
			if ( (h[k-1] < h[k]) && (h[k+1] < h[k]))
			{
				++modes;
				if(modes > PEAK_NUM)
				{
					return false;
				}
			}
		}

		return (modes == PEAK_NUM);
	}

	/*
	* Compute the Intermodes's threshold
	*/
	template<class THistogram, class TOutput>
	void
		IntermodesThresholdCalculator<THistogram, TOutput>
		::GenerateData(void)
	{
		const HistogramType * histogram = this->GetInput();
		// histogram->Print(std::cout);
		if ( histogram->GetTotalFrequency() == 0 )
		{
			itkExceptionMacro(<< "Histogram is empty");
		}
		SizeValueType size = histogram->GetSize(0);

		ProgressReporter progress(this, 0, size );
		if( size == 1 )
		{
			this->GetOutput()->Set( static_cast<OutputType>( histogram->GetMeasurement(0,0) ) );
			return;
		}

		// smooth the histogram
		std::ofstream ofs("IntermodesHistogram.txt");
		std::vector<double> smoothedHist(size);
		for( InstanceIdentifier i = 0; i<size; i++)
		{ 
			smoothedHist[i] = static_cast< double >( histogram->GetFrequency(i, 0) );
			ofs<<  smoothedHist[i]<<"\t";
			progress.CompletedPixel();
		}

		smoothedHist[0] = 0;  // consider many 0 exist.
		ofs<< std::endl;

		SizeValueType smIter = 0;

		while (!BimodalTest(smoothedHist))
		{
			// smooth with a 3 point running mean
			double previous = 0.;
			double current = 0.;
			double next = smoothedHist[0];

			for (size_t i = 0; i < smoothedHist.size() - 1; i++)
			{
				previous = current;
				current = next;
				next = smoothedHist[i + 1];
				smoothedHist[i] = (previous + current + next) / 3.;
			}
			smoothedHist[smoothedHist.size() - 1] = (current + next) / 3.;
			++smIter;

			if (smIter > m_MaximumSmoothingIterations )
			{
				itkGenericExceptionMacro( << "Exceeded maximum iterations for histogram smoothing." );
				return;
			}
		}
		std::cout<< "Converged after "<< smIter << "iterations." <<std::endl;

		double sum = 0;
		for( InstanceIdentifier i = 0; i<size; i++)
		{ 
			sum += smoothedHist[i];
		}
		
		for( InstanceIdentifier i = 0; i<size; i++)
		{ 
			ofs<<  smoothedHist[i] / sum<<"\t";
		}
		ofs<< std::endl;
		for( InstanceIdentifier i = 0; i<size; i++)
		{ 
			ofs<<  histogram->GetMeasurement( i, 0 )<<"\t";
		}
		ofs<< std::endl;
		ofs.close();

		size_t tt = 0;

		if (m_UseInterMode)
		{
			// The threshold is the mean between the two peaks.
			for (size_t i=1; i<smoothedHist.size() - 1; i++)
			{
				if ( ( smoothedHist[i-1] < smoothedHist[i] )&& ( smoothedHist[i+1] < smoothedHist[i] ) )
				{
					tt += i;
				}
			}
			tt /= 2;
		}
		else
		{
			size_t firstpeak = 0;
			size_t secondpeak = 0;
			size_t lowpeak = 0;
			double peak_val1 = 0;
			double  peak_val2 = 0;
			double low_val = 0;

			for (size_t i=1; i<smoothedHist.size() - 1; i++)
			{
				if( (smoothedHist[i-1] < smoothedHist[i] ) && ( smoothedHist[i+1] < smoothedHist[i] ) )
				{
					if( firstpeak == 0)
					{
						firstpeak = i;
						peak_val1 = smoothedHist[i];
					}
					else
					{
						secondpeak = i;
						peak_val2 = smoothedHist[i];
						break;
					}
				}
				
				if( (smoothedHist[i-1] > smoothedHist[i] ) && ( smoothedHist[i+1] > smoothedHist[i] ) )
				{
					lowpeak = i;
					low_val = smoothedHist[i];
				} 
			}

			std::cout<< firstpeak <<"\t"<< lowpeak<< "\t"<< secondpeak<<std::endl;

			double peakRatio = peak_val1 / peak_val2;
			std::cout<<  "Peak ratio"<< peakRatio<<std::endl;

			if(  peakRatio < 100)
			{
				tt = lowpeak;
				std::cout<< "Two peaks: "<< tt<<std::endl;
			}
			else
			{
				for( size_t i = firstpeak; i < smoothedHist.size() ; i++)
				{
					if( smoothedHist[i] < 3 * peak_val2)
					{
						tt = i;
						break;
					}
				}
				std::cout<< "One peak: "<< tt<<std::endl;

				//double total = 0;
				//for ( size_t i=0; i<smoothedHist.size() ; i++)
				//{
				//	total += smoothedHist[i];
				//}
				//double tmpTotal = 0;
				//for (size_t i=smoothedHist.size()  - 1; i >= 0; --i)
				//{
				//	if(  tmpTotal <= total * 0.003)
				//	{
				//		tmpTotal += smoothedHist[i];
				//		tt = i;
				//	}
				//	else
				//	{
				//		std::cout<< tmpTotal <<"\t" << total * 0.005<<std::endl;
				//		std::cout<< "tt2: "<< tt <<std::endl;
				//		break;
				//	}
				//}
			}
		}
		this->GetOutput()->Set( static_cast<OutputType>( histogram->GetMeasurement( tt, 0 ) ) );
	}

	template<class THistogram, class TOutput>
	void
		IntermodesThresholdCalculator<THistogram, TOutput>
		::PrintSelf( std::ostream& os, Indent indent ) const
	{
		Superclass::PrintSelf(os,indent);

		os << indent << "MaximumSmoothingIterations: " << m_MaximumSmoothingIterations << std::endl;
		os << indent << "UseInterMode: " << m_UseInterMode << std::endl;
	}

} // end namespace itk

#endif
