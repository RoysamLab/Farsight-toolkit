//#include "itkRescaleIntensityImageFilter.h"
//#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkNeighborhoodAlgorithm.h"
//---------------------------------------------------------------------------------------------------
//CoOccurrenceFilter Template Begins: This filter is designed to calculate co-occurrence matrix, which
// returns a m_CoOccurrence of size 256 * 256 * 13 
//---------------------------------------------------------------------------------------------------
namespace itk{
template<class TImageType>
class ITK_EXPORT CoOccurrenceFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef CoOccurrenceFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro(CoOccurrenceFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;
  typedef typename Image<unsigned int, 3> COImageType; 
  typedef typename COImageType::IndexType COIndexType;
  typedef typename TImageType::PixelType PixelType;

  typedef typename NumericTraits< PixelType >::AccumulateType SumType;
  typedef typename NumericTraits< SumType >::RealType MeanType;
  typedef typename TImageType::IndexType IndexType;
  typedef typename TImageType::RegionType RegionType;
  
  typedef typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TImageType > FaceCalculatorType;
  typedef typename COImageType::Pointer COType;

  itkGetMacro(CoOccurrence,COType);
  itkSetMacro(CoOccurrence,COType);
  
protected:

	  CoOccurrenceFilter();

protected:

  typedef ConstNeighborhoodIterator<TImageType> NeighborhoodIteratorType;
  typedef ImageRegionIterator<COImageType> COIteratorType;
  void GenerateData();

private:
  CoOccurrenceFilter(Self&);
  void operator=(const Self&);

  COType             m_CoOccurrence;
 };
}

namespace itk
{
 template <class TImageType>
 CoOccurrenceFilter<TImageType>  
 ::CoOccurrenceFilter()
 {

 }

template <class TImageType>
void
CoOccurrenceFilter<TImageType>::
GenerateData()
{ 
  COImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  COImageType::SizeType size;
  size[0]  = 256;
  size[1]  = 256;
  size[2]  = 13;

  COImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);
  
  m_CoOccurrence = COImageType::New();
  m_CoOccurrence->SetRegions(region);
  m_CoOccurrence->Allocate();

  unsigned int element_radius = 3;
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(element_radius);

  FaceCalculatorType faceCalculator;
  FaceCalculatorType::FaceListType faceList;
  faceList = faceCalculator(this->GetInput(), this->GetInput()->GetRequestedRegion(),
  radius);
  FaceCalculatorType::FaceListType::iterator fit;
   // optimizing iteration speed
  for( fit = faceList.begin(); fit != faceList.end(); ++fit)
  {
  NeighborhoodIteratorType it = NeighborhoodIteratorType(radius, this->GetInput(), *fit);

  // NeighborhoodIteratorType it( radius, this->GetInput(), this->GetInput()->GetRequestedRegion());
  NeighborhoodIteratorType::OffsetType offset1;  // offset used for 13 different 3D directions
  NeighborhoodIteratorType::OffsetType offset2;  // offset used for 13 different 3D directions
  
  // initialize Co-Occurrence matrix
  COIteratorType out(m_CoOccurrence,m_CoOccurrence->GetRequestedRegion());
  for(out.GoToBegin(); !out.IsAtEnd(); ++out)
   {
  out.Set(0);
   }
 
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
   { 
   for(int i = 1; i<=13; i++)
	{
	 if(i == 1)
	 {
	 offset1[0] = 1;
	 offset1[1] = 0;
	 offset1[2] = 0;
	 offset2[0] = -1;
	 offset2[1] = 0;
	 offset2[2] = 0;
	 }
	 if(i == 2)
	 {
     offset1[0] = 1;
	 offset1[1] = 0;
	 offset1[2] = 1;
	 offset2[0] = -1;
	 offset2[1] = 0;
	 offset2[2] = -1;
	 }
	 if(i == 3)
	 {
	 offset1[0] = 0;
	 offset1[1] = 0;
	 offset1[2] = -1;
	 offset2[0] = 0;
	 offset2[1] = 0;
	 offset2[2] = 1;
	 }
	 if(i == 4)
	 {
	 offset1[0] = -1;
	 offset1[1] = 0;
	 offset1[2] = 1;
	 offset2[0] = 1;
	 offset2[1] = 0;
	 offset2[2] = -1;
	 }
	 if(i == 5)
	 {
	 offset1[0] = 0;
	 offset1[1] = -1;
	 offset1[2] = 1;
	 offset2[0] = 0;
	 offset2[1] = 1;
	 offset2[2] = -1;
	 }
	 if(i == 6)
	 {
	 offset1[0] = 0;
	 offset1[1] = -1;
	 offset1[2] = 0;
	 offset2[0] = 0;
	 offset2[1] = 1;
	 offset2[2] = 0;
	 }
	 if(i == 7)
	 {
	 offset1[0] = 0;
	 offset1[1] = 1;
	 offset1[2] = 1;
	 offset2[0] = 0;
	 offset2[1] = -1;
	 offset2[2] = -1;
	 }
	 if(i == 8)
	 {
	 offset1[0] = 1;
	 offset1[1] = 1;
	 offset1[2] = 0;
	 offset2[0] = -1;
	 offset2[1] = -1;
	 offset2[2] = 0;
	 }
	 if(i == 9)
	 {
	 offset1[0] = 1;
	 offset1[1] = -1;
	 offset1[2] = 0;
	 offset2[0] = -1;
	 offset2[1] = 1;
	 offset2[2] = 0;
	 }
	 if(i == 10)
	 {
	 offset1[0] = -1;
	 offset1[1] = 1;
	 offset1[2] = -1;
	 offset2[0] = 1;
	 offset2[1] = -1;
	 offset2[2] = 1;
	 }
	 if(i == 11)
	 {
	 offset1[0] = -1;
	 offset1[1] = -1;
	 offset1[2] = -1;
	 offset2[0] = 1;
	 offset2[1] = 1;
	 offset2[2] = 1;
	 }
	 if(i == 12)
	 {
	 offset1[0] = 1;
	 offset1[1] = 1;
	 offset1[2] = -1;
	 offset2[0] = -1;
	 offset2[1] = -1;
	 offset2[2] = 1;
	 }
	 if(i == 13)
	 {
	 offset1[0] = 1;
	 offset1[1] = -1;
	 offset1[2] = -1;
	 offset2[0] = -1;
	 offset2[1] = 1;
	 offset2[2] = 1;
	 }

     COIndexType index1;
     index1[0]  = static_cast<unsigned int>(it.GetCenterPixel());
	 index1[1]  = static_cast<unsigned int>(it.GetPixel(offset1));
	 index1[2]  = i - 1;
     COIndexType index2;
	 index2[0]  = static_cast<unsigned int>(it.GetCenterPixel());
	 index2[1]  = static_cast<unsigned int>(it.GetPixel(offset2));
	 index2[2]  = i - 1;

     m_CoOccurrence->SetPixel(index1, static_cast<unsigned int>(m_CoOccurrence->GetPixel(index1) + 1));  // m_CoOccurrence is 256 * 256 * 13 image containing
	 m_CoOccurrence->SetPixel(index2, static_cast<unsigned int>(m_CoOccurrence->GetPixel(index2) + 1));  // co-occurrence matrixs for 13 directions
	}
   }
  }
  	 std::cout<<"Obtained Co-Occurrence Matrixs for 13 directions"<<std::endl;
  
 
  
  this->GraftOutput( this->GetOutput() );  
  this->GraftOutput(this->GetOutput());
    
}


template <class TImageType>
void
CoOccurrenceFilter<TImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "m_CoOccurrence" << this->m_CoOccurrence
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//CoOccurrenceFilter Template Ends
//---------------------------------------------------------------------------------------------------