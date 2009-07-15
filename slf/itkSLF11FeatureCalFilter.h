///////////////////////////////////////////////////////////////////////////////////////////////////////////
////         SLF11 Feature Calculation Filter, Translation of Murphy's SLF code                       /////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
 
#include "itkImageToImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
      
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkVector.h"


//---------------------------------------------------------------------------------------------------
//SLF11FeatureCalFilter Template
//---------------------------------------------------------------------------------------------------
namespace itk{
template<class TImageType, class TImageType1>
class ITK_EXPORT SLF11FeatureCalFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef SLF11FeatureCalFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro(SLF11FeatureCalFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;
   
  typedef typename TImageType::PixelType PixelType;
  
  typedef typename NumericTraits< PixelType >::AccumulateType SumType;
  typedef typename NumericTraits< SumType >::RealType MeanType;
  typedef typename TImageType::IndexType IndexType;
  typedef typename Vector<float, 28> VectorType;

  itkGetMacro(FeatureVector,VectorType); 
  itkSetMacro(FeatureVector,VectorType);

  itkGetMacro(Para, int);  
  itkSetMacro(Para, int);             
  
  itkGetMacro(Variance, float);
  itkSetMacro(Variance, float);

protected:

	  SLF11FeatureCalFilter();

protected:

  typedef ImageRegionConstIterator<TImageType> ConstIteratorType;
  typedef ImageRegionIterator<TImageType> IteratorType;
  typedef ImageRegionConstIteratorWithIndex<TImageType> ConstIteratorWithIndexType;
   
  typedef EdgeFeaturesFilter<TImageType,TImageType1> EdgeFeaturesFilterType;
  typedef TextureFeaturesFilter<TImageType> TextureFeaturesFilterType;
  void GenerateData();

private:
  SLF11FeatureCalFilter(Self&);
  void operator=(const Self&);
  
  int                         m_Para;
  VectorType         m_FeatureVector;
  float                   m_Variance;
 };
}

namespace itk
{
 template <class TImageType, class TImageType1>
 SLF11FeatureCalFilter<TImageType,TImageType1>  
 ::SLF11FeatureCalFilter()
 {
 }

template <class TImageType, class TImageType1>
void
SLF11FeatureCalFilter<TImageType,TImageType1>::
GenerateData()
{     
   m_Para = 1;
   m_Variance = 1;

   EdgeFeaturesFilterType::Pointer filter1 = EdgeFeaturesFilterType::New();
   TextureFeaturesFilterType::Pointer filter2 = TextureFeaturesFilterType::New();
   filter1->SetInput( this->GetInput() );
   filter1->SetVariance(m_Variance);
//  filter1->SetMaskInput( filter->GetOutput() );
  filter1->Update();
  std::cout<<"Finished calculating Edge Features"<<std::endl;
  filter2->SetInput( this->GetInput() );

  m_FeatureVector[0] = filter1->GetPenpixel();
  m_FeatureVector[1] = filter1->GetPenflu();
  filter2->Update();
  for(int i = 0; i < 26; i++)
  {
  m_FeatureVector[ i + 2] = filter2->GetFeatureVector()[i];
  }

  filter1->GraftOutput( this->GetOutput() );
  this->GraftOutput(filter1->GetOutput());
    
}


template <class TImageType, class TImageType1>
void
SLF11FeatureCalFilter<TImageType,TImageType1>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "Threshold:" <<m_Para
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//SLF11FeatureCalFilter Template
//---------------------------------------------------------------------------------------------------