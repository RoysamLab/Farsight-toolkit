///////////////////////////////////////////////////////////////////////////////////////////////////////////
////         Finding Objects in the 3D image, Translation of Murphy's SLF code "ml_3dfindobj.c"      /////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "itkRescaleIntensityImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"

#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkVector.h"

#include "itkCastImageFilter.h"
//---------------------------------------------------------------------------------------------------
//EdgeFeaturesFilter Template Begins
//---------------------------------------------------------------------------------------------------
namespace itk{
template<class TImageType, class TImageType1>
class ITK_EXPORT EdgeFeaturesFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef EdgeFeaturesFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro(EdgeFeaturesFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;
   
  typedef typename TImageType::PixelType PixelType;
  typedef typename Vector<float, 2> VectorType;

  typedef typename NumericTraits< PixelType >::AccumulateType SumType;
  typedef typename NumericTraits< SumType >::RealType MeanType;
  typedef typename TImageType::IndexType IndexType;
  typedef typename TImageType::RegionType RegionType;
  
 // itkGetMacro(FeatureVector,VectorType); 
 // itkSetMacro(FeatureVector,VectorType);

  itkGetMacro(Penpixel,MeanType);
  itkSetMacro(Penpixel,MeanType);

  itkGetMacro(Penflu,MeanType);
  itkSetMacro(Penflu,MeanType);

  itkGetMacro(Variance, float);
  itkSetMacro(Variance, float);
  
  
protected:

	  EdgeFeaturesFilter();

protected:

  typedef ImageRegionConstIterator<TImageType> ConstIteratorType;

  typedef ImageRegionIterator<TImageType> IteratorType;

  typedef CannyEdgeDetectionImageFilter<TImageType1,TImageType1> CannyFilterType;
  
  typedef RescaleIntensityImageFilter< TImageType, TImageType > RescalerType;
  
  typedef BinarizeFilter<TImageType> BinarizeFilterType;
  
  typedef CastImageFilter<TImageType, TImageType1> CasterType;
  typedef CastImageFilter<TImageType1, TImageType> CasterType1;

  void GenerateData();

private:
  EdgeFeaturesFilter(Self&);
  void operator=(const Self&);
  typename RescalerType::Pointer      m_RescaleFilter;
  typename CasterType::Pointer        m_CasterFilter;
  typename CasterType1::Pointer       m_CasterFilter1;

  VectorType             m_FeatureVector;
  MeanType               m_Penflu;
  MeanType               m_Penpixel;
  float                  m_Variance;


 };
}

namespace itk
{
 template <class TImageType, class TImageType1>
 EdgeFeaturesFilter<TImageType, TImageType1>  
 ::EdgeFeaturesFilter()
 {
	 m_CasterFilter = CasterType::New();
	 m_CasterFilter1 = CasterType1::New();
 }

template <class TImageType, class TImageType1>
void
EdgeFeaturesFilter<TImageType, TImageType1>::
GenerateData()
{     
  m_Variance = 1;
 
  CannyFilterType::Pointer filter1 = CannyFilterType::New();                       //Segmentation Filter
  BinarizeFilterType::Pointer filter2 = BinarizeFilterType::New();
  
 
  filter2->SetInput(this->GetInput());
 // filter2->Update();

  m_CasterFilter->SetInput(filter2->GetOutput());
  m_CasterFilter->Update();

  filter1->SetInput(m_CasterFilter->GetOutput());
  filter1->SetVariance(m_Variance);
  m_CasterFilter1->SetInput(filter1->GetOutput());
 // filter1->Update();
  m_CasterFilter1->Update();
  	ConstIteratorType maskIt( filter2->GetOutput(), filter2->GetOutput()->GetRequestedRegion());
	ConstIteratorType fluIt(  this->GetInput(), this->GetInput()->GetRequestedRegion());
    ConstIteratorType edgeIt(  m_CasterFilter1->GetOutput(), m_CasterFilter1->GetOutput()->GetRequestedRegion());
    
	SumType pixelt = NumericTraits< SumType >::Zero;
	SumType fluort = NumericTraits< SumType >::Zero;
    m_Penpixel = NumericTraits< MeanType >::Zero;
    m_Penflu = NumericTraits< MeanType >::Zero;

    for ( maskIt.GoToBegin(),fluIt.GoToBegin(); !maskIt.IsAtEnd(),!fluIt.IsAtEnd(); ++maskIt, ++fluIt)
    {
	  if(maskIt.Get() == 1) 
	  {
	 // maskIt.Set(1);
	  pixelt += 1;
	  fluort += fluIt.Get();
	  }
    }

	for ( maskIt.GoToBegin(), fluIt.GoToBegin(), edgeIt.GoToBegin(); !maskIt.IsAtEnd(), !fluIt.IsAtEnd(),!edgeIt.IsAtEnd(); ++maskIt,++fluIt, ++edgeIt)
    {
      if(edgeIt.Get() == 1 && maskIt.Get() ==1 )
	  {
	//  edgeIt.Set(1); 
      m_Penpixel += 1;
	  m_Penflu += fluIt.Get();
	  }
    }

  m_Penpixel = m_Penpixel/pixelt * 100;
  m_Penflu   = m_Penflu/fluort * 100;

//  m_FeatureVector[0] = m_Penpixel;
//  m_FeatureVector[1] = m_Penflu;

  filter2->GraftOutput( this->GetOutput() );  
  this->GraftOutput(filter2->GetOutput());
    
}


template <class TImageType, class TImageType1>
void
EdgeFeaturesFilter<TImageType, TImageType1>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "Threshold:" << this->m_Variance
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//EdgeFeaturesFilter Template Ends
//---------------------------------------------------------------------------------------------------