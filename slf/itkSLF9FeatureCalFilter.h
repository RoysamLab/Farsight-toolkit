///////////////////////////////////////////////////////////////////////////////////////////////////////////
////         SLF9 Feature Calculation Filter, Translation of Murphy's SLF code                       /////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Input: Threshold
// Output: FeatureVector containing all 28 features
// Use: 
//   Binarization Filter: itkBinarizeFilter.h  
//   Objects Segmentation Filter: itkFindObjFilter.h
//   Protein Feature Calculation Filter: itkObj2FeatureFilter.h
//   DNA Feature Calculation Filter: itkObj2FeatureDNAFilter.h

 
#include "itkImageToImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
      
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkVector.h"


//---------------------------------------------------------------------------------------------------
//SLF9FeatureCalFilter Template
//---------------------------------------------------------------------------------------------------
namespace itk{
template<class TImageType, class TDNAImage>
class ITK_EXPORT SLF9FeatureCalFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef SLF9FeatureCalFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro(SLF9FeatureCalFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;
   
  typedef typename TImageType::PixelType PixelType;
  
  typedef typename NumericTraits< PixelType >::AccumulateType SumType;
  typedef typename NumericTraits< SumType >::RealType MeanType;
  typedef typename TImageType::IndexType IndexType;
  typedef typename Vector<float, 28> VectorType;
  typedef typename Vector<float, 3> COFVectorType;
  itkGetMacro(FeatureVector,VectorType); 
  itkSetMacro(FeatureVector,VectorType);

  itkGetMacro(Threshold, int);  
  itkSetMacro(Threshold, int);             
  
  void SetDNAInput(TDNAImage *input)
    {
      this->SetNthInput(1, const_cast<TDNAImage *>(input) );
    }

  TDNAImage * GetDNAInput()
    {
      return static_cast<TDNAImage*>(const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
    }

protected:

	  SLF9FeatureCalFilter();

protected:

  typedef ImageRegionConstIterator<TImageType> ConstIteratorType;
  typedef ImageRegionIterator<TImageType> IteratorType;
  typedef ImageRegionConstIteratorWithIndex<TImageType> ConstIteratorWithIndexType;
  
  typedef BinarizeFilter<TImageType> BinarizeFilterType;
  typedef FindObjFilter<TImageType> FindObjFilterType;
  typedef Obj2FeatureFilter<TImageType> Obj2FeatureFilterType;
  typedef Obj2FeatureDNAFilter<TImageType> Obj2FeatureDNAFilterType;
  void GenerateData();

private:
  SLF9FeatureCalFilter(Self&);
  void operator=(const Self&);

  int                    m_Threshold;
  VectorType    m_FeatureVector;
 };
}

namespace itk
{
 template <class TImageType, class TDNAImage>
 SLF9FeatureCalFilter<TImageType, TDNAImage>  
 ::SLF9FeatureCalFilter()
 {
 }

template <class TImageType, class TDNAImage>
void
SLF9FeatureCalFilter<TImageType,TDNAImage>::
GenerateData()
{     
  
  BinarizeFilterType::Pointer filter = BinarizeFilterType::New();               //Binarization Filter
  BinarizeFilterType::Pointer filterDNA = BinarizeFilterType::New();               //Binarization Filter
  FindObjFilterType::Pointer filter1 = FindObjFilterType::New();                //Objects Segmentation Filter
  FindObjFilterType::Pointer filter1DNA = FindObjFilterType::New();                //Objects Segmentation Filter
  Obj2FeatureFilterType::Pointer filter2 = Obj2FeatureFilterType::New();        //Protein Feature Calculation Filter
  Obj2FeatureDNAFilterType::Pointer filter2DNA = Obj2FeatureDNAFilterType::New();  //DNA Feature Calculation Filter
  

// Calculate Protein Image related features
  filter->SetInput( this->GetInput() );
  filter1->SetInput(filter->GetOutput());
  filter1->SetThreshold(m_Threshold);
  filter1->Update();

  COFVectorType COFProt = filter1->GetCOF();
  float     VolProt = filter1->GetAverageVolume() * filter1->GetNobject();
  filter2->SetInput(this->GetInput());
  filter2->SetNobject(filter1->GetNobject());
  filter2->SetAverageVolume(filter1->GetAverageVolume());
  filter2->SetSDVolume(filter1->GetSDVolume());
  filter2->SetMax2Min(filter1->GetMax2Min());
  filter2->SetCOF(filter1->GetCOF());
  filter2->SetBlob(filter1->GetBlob());
  filter2->Update();

  VectorType FeatureVectorProt = filter2->GetFeatureVector();
   //  std::cout<<"----------------------------SLF9 Feature Set----------------------------"<<std::endl;
   //for(int i=0; i<28; i++) 
   //{
   //	 std::cout<<"SLF9_"<<i+1<<":"<<FeatureVectorProt[i]<<std::endl;
   //}
   //std::cout<<"----------------------------SLF9 Feature Set----------------------------"<<std::endl;



 // Calculate DNA Image related features
  filterDNA->SetInput(this->GetDNAInput());
  filter1DNA->SetInput(filterDNA->GetOutput());
  filter1DNA->SetThreshold(m_Threshold);
  filter1DNA->Update();
  
  TImageType::Pointer DNALabelImage = filter1DNA->GetOutput();

  filter2DNA->SetCOFProt(COFProt);
  filter2DNA->SetVolProt(filter1->GetTotalVolume());
  
  filter2DNA->SetInput(this->GetDNAInput());
  filter2DNA->SetNobjectDNA(filter1->GetNobject());       // using objects from Protein images instead of DNA images
  filter2DNA->SetAverageVolumeDNA(filter1DNA->GetAverageVolume());
  filter2DNA->SetTotalVolumeDNA(filter1DNA->GetTotalVolume());
  filter2DNA->SetCOFDNA(filter1DNA->GetCOF());
  filter2DNA->SetBlobDNA(filter1->GetBlob());             // using objects blob from Protein images instead of DNA images
  filter2DNA->Update();


 // fraction of protein fluorescence overlapping with DNA fluorescence 
  ConstIteratorType it( filter->GetOutput(),filter->GetOutput()->GetRequestedRegion());
  ConstIteratorType ot( filterDNA->GetOutput(),filterDNA->GetOutput()->GetRequestedRegion() );
   int fractionNumber = 0;
   for(it.GoToBegin(),ot.GoToBegin(); !it.IsAtEnd(), !ot.IsAtEnd(); ++it, ++ot)
   {
    if(it.Get()!=0 && ot.Get()!=0)
    {
     fractionNumber += 1;
    }
   }
  float fraction = fractionNumber/VolProt;

  VectorType FeatureVectorDNA = filter2DNA->GetFeatureVector();
  m_FeatureVector    = FeatureVectorProt + FeatureVectorDNA;
  m_FeatureVector[13]           = fraction;

    
  // std::cout<<"----------------------------SLF9 Feature Set----------------------------"<<std::endl;
  // for(int i=0; i<28; i++) 
  // {
  //	 std::cout<<"SLF9_"<<i+1<<":"<<m_FeatureVector[i]<<std::endl;
  // }
  // std::cout<<"----------------------------SLF9 Feature Set----------------------------"<<std::endl;

  filterDNA->GraftOutput( this->GetOutput() );
  this->GraftOutput(filterDNA->GetOutput());
    
}


template <class TImageType, class TDNAImage>
void
SLF9FeatureCalFilter<TImageType,TDNAImage>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "Threshold:" << this->m_Threshold
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//SLF9FeatureCalFilter Template
//---------------------------------------------------------------------------------------------------