
#include "itkVector.h"


//---------------------------------------------------------------------------------------------------
//TextureFeaturesFilter Template Begins
//---------------------------------------------------------------------------------------------------
namespace itk{
template<class TImageType>
class ITK_EXPORT TextureFeaturesFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef TextureFeaturesFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro(TextureFeaturesFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;
   
  typedef typename TImageType::PixelType PixelType;
  typedef typename Vector<float, 26> VectorType;

  typedef typename NumericTraits< PixelType >::AccumulateType SumType;
  typedef typename NumericTraits< SumType >::RealType MeanType;
  typedef typename TImageType::IndexType IndexType;
  typedef typename TImageType::RegionType RegionType;
  
 // itkGetMacro(FeatureVector,VectorType); 
 // itkSetMacro(FeatureVector,VectorType);

  itkGetMacro(FeatureVector,VectorType);
  itkSetMacro(FeatureVector,VectorType);
  
protected:

	  TextureFeaturesFilter();

protected:
  typedef CoOccurrenceFilter<TImageType> CoOccurrenceFilterType;
  typedef CO2Feature3DFilter<TImageType> CO2Feature3DFilterType;

  void GenerateData();

private:
  TextureFeaturesFilter(Self&);
  void operator=(const Self&);

  VectorType             m_FeatureVector;

 };
}

namespace itk
{
 template <class TImageType>
 TextureFeaturesFilter<TImageType>  
 ::TextureFeaturesFilter()
 {

 }

template <class TImageType>
void
TextureFeaturesFilter<TImageType>::
GenerateData()
{     
  CoOccurrenceFilterType::Pointer filter = CoOccurrenceFilterType::New();
  CO2Feature3DFilterType::Pointer filter1  = CO2Feature3DFilterType::New();
  filter->SetInput( this->GetInput() );
  filter->Update();

  filter1->SetInput( this->GetInput() );
  filter1->SetCoOccurrence(filter->GetCoOccurrence());  // Get Occurrence matrixs from filter
  filter1->Update();
  m_FeatureVector = filter1->GetFeatureVector();
  
  this->GraftOutput( this->GetOutput() );  
    
}


template <class TImageType>
void
TextureFeaturesFilter<TImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "FeatureVector:" << this->m_FeatureVector
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//TextureFeaturesFilter Template Ends
//---------------------------------------------------------------------------------------------------