///////////////////////////////////////////////////////////////////////////////////////////////////////////
////Background Substraction. Translation of Murphy's Code "ml_binarize.c", using NIH Threshold Method/////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <math.h>

namespace itk{
template<class TImageType>
class ITK_EXPORT BinarizeFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef BinarizeFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(BinarizeFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;

  typedef typename TImageType::PixelType PixelType;

  itkGetMacro(Para, int);
  itkSetMacro(Para, int);

protected:

	  BinarizeFilter();

protected:
  
  typedef itk::ImageRegionConstIterator< TImageType > ConstIteratorType;
  typedef ImageRegionIterator< TImageType> IteratorType;
  typedef RescaleIntensityImageFilter< TImageType, TImageType > RescalerType;
  void GenerateData();

private:
  BinarizeFilter(Self&);
  void operator=(const Self&);
  
  typename RescalerType::Pointer     m_RescaleFilter1;
   int m_Para;


 };
}

namespace itk
{
 template <class TImageType>
 BinarizeFilter<TImageType>  
 ::BinarizeFilter()
 {
     m_RescaleFilter1 = RescalerType::New();	 
	 m_Para = 1;

	 m_RescaleFilter1->SetOutputMinimum(0);
     m_RescaleFilter1->SetOutputMaximum(255);
//	 std::cout<<NumericTraits<PixelType>::max()<<NumericTraits<PixelType>::NonpositiveMin()<<std::endl;

 }

template <class TImageType>
void
BinarizeFilter<TImageType>::
GenerateData()
{ 
  unsigned long *hist = new unsigned long [256];
  for(int i=0; i<256; i++){
  hist[i] = 0;
  }
  int Threshold = 0;
  int Result    = 0;

  TImageType::Pointer output = TImageType::New();
  TImageType::RegionType inputRegion  = this->GetInput()->GetRequestedRegion();
  TImageType::RegionType outputRegion = this->GetInput()->GetRequestedRegion();

  output->SetRegions(outputRegion);
  output->Allocate();
  ConstIteratorType inputIt(  this->GetInput(), inputRegion  );
  IteratorType      outputIt(  output,         outputRegion );
  
  for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
    { 
     hist[ static_cast<int>(inputIt.Get())] += 1;
    }
//  for(int k=0; k<256; k++)
// {
//	  std::cout<<"Histogram:"<<hist[k]<<std::endl;
// }
  
  
// NIH Threshold Method
  hist[0]   = 0;
  hist[255] = 0;
  int error = 0;
  int MinIndex = 0;
  int MaxIndex = 255;

  while(hist[MinIndex]==0 && MinIndex < 255){
  MinIndex += 1; 
  }
  while(hist[MaxIndex]==0 && MaxIndex > 0){
  MaxIndex -= 1; 
  }

  if(MinIndex >= MaxIndex){
  Result = 128;
  std::cout<<"There is # <=1 nonzero pixel intensity "<<std::endl;
  error  = 1; 
  }

  int MovingIndex = MinIndex;
  if(error == 0){
     Result = MaxIndex;
	// std::cout<<"Result=MaxIndex:"<<Result<<std::endl;
	 while((MovingIndex + 1) <= Result){
	   if(MovingIndex> (MaxIndex-1)){
		   break;
	   }
	   unsigned long sum1 = 0;
       unsigned long sum2 = 0;
       unsigned long sum3 = 0;
       unsigned long sum4 = 0;

	   for(int i = MinIndex; i <= MovingIndex; i++){
	   sum1 += i * hist[i];
	   }
	   for(int i = MinIndex; i <= MovingIndex; i++){
	   sum2 += hist[i];
	   }
	   for(int i = MovingIndex + 1; i <= MaxIndex; i++){
	   sum3 += i * hist[i];
	   }
	   for(int i = MovingIndex + 1; i <= MaxIndex; i++){
	   sum4 += hist[i];
	   }
	   
	   Result = (sum1/sum2 + sum3/sum4)/2;
	   MovingIndex = MovingIndex + 1;
	 }
  }
 delete[] hist;
	  
//  Threshold = 255 - static_cast<unsigned int>(Result) + 1;
    Threshold =  Result;

	std::cout<<"Binarization Threshold:"<<Threshold<<std::endl;
  for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
        ++inputIt, ++outputIt)
    {
	int newvalue = inputIt.Get() ;
	if ( newvalue < Threshold ) newvalue = 0;
	else newvalue = 1;
    outputIt.Set(newvalue);
    }

//  output->CopyInformation(m_RescaleFilter1->GetOutput());
  
 // m_RescaleFilter2->SetInput(output);
 
  //output->GraftOutput( this->GetOutput() );
 // m_RescaleFilter2->Update();
   
  this->GraftOutput(output);
    
}


template <class TImageType>
void
BinarizeFilter<TImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "Element Radius:" << this->m_Para
    << std::endl;
}

}
