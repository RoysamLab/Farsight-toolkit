///////////////////////////////////////////////////////////////////////////////////////////////////////////
////Background Substraction. Translation of Murphy's Code "ml_3dbgsub.c", using NIH Threshold Method/////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
//#include <math.h>

//---------------------------------------------------------------------------------------------------
//TDBGSUBFilter Template Begins
//---------------------------------------------------------------------------------------------------

namespace itk{
template<class TImageType>
class ITK_EXPORT TDBGSUBFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef TDBGSUBFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(TDBGSUBFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;

  typedef typename TImageType::PixelType PixelType;
//  itkGetMacro( Threshold, PixelType);
// itkSetMacro( Threshold, PixelType);
  itkGetMacro(Para, int);
  itkSetMacro(Para, int);

protected:

	  TDBGSUBFilter();

protected:
  
  typedef itk::ImageRegionConstIterator< TImageType > ConstIteratorType;
  typedef ImageRegionIterator< TImageType> IteratorType;
  typedef RescaleIntensityImageFilter< TImageType, TImageType > RescalerType;
  void GenerateData();

private:
  TDBGSUBFilter(Self&);
  void operator=(const Self&);
  
 // typename ShapedNeighborhoodIteratorType::Pointer m_ShapedNeighborhoodIterator;
 // typename IteratorType::Pointer                   m_Iterator;

  typename RescalerType::Pointer     m_RescaleFilter2;
   int m_Para;
//  TImageType inputImage;

 };
}

namespace itk
{
 template <class TImageType>
 TDBGSUBFilter<TImageType>  
 ::TDBGSUBFilter()
 {
	// ShapedNeighborhoodIteratorType::Pointer m_ShapedNeighborhoodIterator = ShapedNeighborhoodIteratorType::New();
	// m_Iterator = IteratorType::New();

	 m_RescaleFilter2 = RescalerType::New();
//	 FaceCalculatorType m_FaceCalculator = FaceCalculatorType::New();   	 
	 m_Para = 1;

//	 std::cout<<NumericTraits<PixelType>::max()<<NumericTraits<PixelType>::NonpositiveMin()<<std::endl;
	 m_RescaleFilter2->SetOutputMinimum(0);
     m_RescaleFilter2->SetOutputMaximum(255);
 }

template <class TImageType>
void
TDBGSUBFilter<TImageType>::
GenerateData()
{
  unsigned long hist[256];
//  m_RescaleFilter1->SetInput( this->GetInput() );
//  m_RescaleFilter1->Update();
///////////////////////////////////////////
  TImageType::Pointer output = TImageType::New();
  TImageType::RegionType inputRegion  = this->GetOutput()->GetRequestedRegion();
  TImageType::RegionType outputRegion = this->GetOutput()->GetRequestedRegion();

  output->SetRegions(outputRegion);
  output->Allocate();

  ConstIteratorType inputIt(  this->GetOutput(), inputRegion  );
  IteratorType      outputIt(  output,         outputRegion );
  
  for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
    { 
     hist[inputIt.Get()]++;
    }

//  unsigned long highest_count = hist[1]; // pixel value in second largest number
  unsigned long highest_count = hist[0];
  int most_common_idx = 0;
  for(int i = 1; i < 256; i++) {
    if( hist[i] > highest_count ) {
      highest_count = hist[i];
      most_common_idx = i;
     }
   }
  std::cout<<hist<<std::endl;
  unsigned long MostCommonPixelValue = most_common_idx;

  for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
        ++inputIt, ++outputIt)
    {
	int newvalue = inputIt.Get() - MostCommonPixelValue;
	if ( newvalue < 0 ) newvalue = 0;
    outputIt.Set(newvalue);
    }
//  output->CopyInformation(m_RescaleFilter1->GetOutput());
  
  m_RescaleFilter2->SetInput(output);
 
//////////////////////////////////////////
  m_RescaleFilter2->GraftOutput( this->GetOutput() );
  m_RescaleFilter2->Update();
   
  this->GraftOutput(m_RescaleFilter2->GetOutput());
  
}


template <class TImageType>
void
TDBGSUBFilter<TImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "Reserved Para:" << this->m_Para
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//TDBGSUBFilter Template Ends
//---------------------------------------------------------------------------------------------------