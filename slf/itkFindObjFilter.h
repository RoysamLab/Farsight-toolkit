///////////////////////////////////////////////////////////////////////////////////////////////////////////
////         Finding Objects in the 3D image, Translation of Murphy's SLF code "ml_3dfindobj.c"      /////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "itkRescaleIntensityImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkBlobSpatialObject.h"
#include "itkSpatialObjectPoint.h"
//#include "itkSceneSpatialObject.h"
#include "itkNumericTraits.h"
#include "itkLabelStatisticsImageFilter.h"
#include <math.h> 


//---------------------------------------------------------------------------------------------------
//FindObjFilter Template Begins
//---------------------------------------------------------------------------------------------------
namespace itk{
template<class TImageType>
class ITK_EXPORT FindObjFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef FindObjFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro(FindObjFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;
   
  typedef typename TImageType::PixelType PixelType;

  typedef BlobSpatialObject<3> BlobType;
  typedef BlobType::Pointer BlobPointer;
  typedef SpatialObjectPoint<3> BlobPointType;
  
  typedef typename Vector<float, 3> COFVectorType;
  typedef typename NumericTraits< PixelType >::AccumulateType SumType;
  typedef typename NumericTraits< SumType >::RealType MeanType;
  typedef typename TImageType::IndexType IndexType;
  typedef typename TImageType::RegionType RegionType;
  

  itkGetMacro(Blob,BlobPointer);
  itkSetMacro(Blob,BlobPointer);

  itkGetMacro(Threshold, int);
  itkSetMacro(Threshold, int);
  
  itkGetMacro(TotalVolume,MeanType);
  itkSetMacro(TotalVolume,MeanType);

  itkGetMacro(AverageVolume,MeanType);
  itkSetMacro(AverageVolume,MeanType);
  itkGetMacro(SDVolume,MeanType);
  itkSetMacro(SDVolume,MeanType);
  itkGetMacro(Max2Min,MeanType);
  itkSetMacro(Max2Min,MeanType);
  
  itkGetMacro(COF, COFVectorType);
  itkSetMacro(COF, COFVectorType);
  itkGetMacro(Nobject, int);
  itkSetMacro(Nobject, int);


  
protected:

	  FindObjFilter();

protected:

  typedef ImageRegionConstIteratorWithIndex<TImageType> ConstIteratorWithIndexType;

  typedef ImageRegionConstIterator<TImageType> ConstIterator;

  typedef ConnectedComponentImageFilter< TImageType,
                           TImageType > ConnectedFilterType;
  typedef RelabelComponentImageFilter< TImageType,
                           TImageType > RelabelFilterType;
  typedef LabelStatisticsImageFilter<TImageType, 
	                       TImageType> LabelStatisticsImageFilterType;
  
  typedef RescaleIntensityImageFilter< TImageType, TImageType > RescalerType;
  
  void GenerateData();

private:
  FindObjFilter(Self&);
  void operator=(const Self&);
  typename RescalerType::Pointer     m_RescaleFilter;

  int                    m_Threshold;
  MeanType               m_AverageVolume;
  MeanType               m_SDVolume;
  float                  m_Max2Min;
  typename BlobPointer   m_Blob;    
  COFVectorType          m_COF;
  int                    m_Nobject;
  MeanType               m_TotalVolume;

 };
}

namespace itk
{
 template <class TImageType>
 FindObjFilter<TImageType>  
 ::FindObjFilter()
 {
 }

template <class TImageType>
void
FindObjFilter<TImageType>::
GenerateData()
{     
     m_TotalVolume =  NumericTraits< MeanType >::Zero;
	 m_Blob    = BlobType::New();
	ConnectedFilterType::Pointer filter1 = ConnectedFilterType::New();                       //Segmentation Filter
    RelabelFilterType::Pointer filter2 = RelabelFilterType::New();                           //Relabel Filter
	LabelStatisticsImageFilterType::Pointer filter3 = LabelStatisticsImageFilterType::New(); //Filter used to query regions by label
    
	filter1->SetInput(this->GetInput());
    filter2->SetInput(filter1->GetOutput());
    filter2->Update();
    
	filter3->SetInput(this->GetInput());
	filter3->SetLabelInput(filter2->GetOutput());
	filter3->Update();

	m_Nobject = filter2->GetNumberOfObjects();
	

	for(int j = 1; j<=m_Nobject; j--)
	{
	m_TotalVolume += filter2->GetSizeOfObjectInPixels(j);
	}
	//std::cout<<"m_TotalVolume"<<m_TotalVolume<<std::endl;
	int i = m_Nobject;
    while(filter2->GetSizeOfObjectInPixels(i) < m_Threshold)
	{
	m_Nobject--;
	i--; 
	}

	std::cout<<"Input Threshold:"<<m_Threshold<<std::endl;
	std::cout<<"Number Of Objects Before Thresholding:"<<filter2->GetNumberOfObjects()<<std::endl;  //No. of objects
    std::cout<<"Number Of Objects After Thresholding:"<<m_Nobject<<std::endl;  //No. of objects
    
	// SLF Feature: (3) Average object volume (average number of above threshold voxels per object)
	m_AverageVolume = NumericTraits< MeanType >::Zero;
	double sum     = 0;
	for (int j = 1; j <= m_Nobject; j++)
	{
	m_AverageVolume += filter2->GetSizeOfObjectInPixels(j)/m_Nobject;
	}

	for (int k = 1; k<= m_Nobject; k++)
	{
	sum += pow((filter2->GetSizeOfObjectInPixels(k) - m_AverageVolume), 2);
	}
    
    m_SDVolume = sqrt(sum / m_Nobject);
	m_Max2Min  = filter2->GetSizeOfObjectInPixels(1)/filter2->GetSizeOfObjectInPixels(m_Nobject);

    
    BlobType::PointListType list;

     SumType TotalSum = NumericTraits< MeanType >::Zero;
     MeanType TotalWeightX = NumericTraits< MeanType >::Zero;
     MeanType TotalWeightY = NumericTraits< MeanType >::Zero;
     MeanType TotalWeightZ = NumericTraits< MeanType >::Zero;
     

	for(int j=1; j<=m_Nobject; j++)
	{
     RegionType labelRegion = filter3->GetRegion(j);
     ConstIteratorWithIndexType it(this->GetInput(), labelRegion);
	 SumType sum1 = NumericTraits< MeanType >::Zero;
     MeanType WeightX = NumericTraits< MeanType >::Zero;
     MeanType WeightY = NumericTraits< MeanType >::Zero;
     MeanType WeightZ = NumericTraits< MeanType >::Zero;
     BlobPointType p;
	 for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	 {
   //   IndexType idx = it.GetIndex();
	  WeightX += it.Get() * it.GetIndex()[0];
	  WeightY += it.Get() * it.GetIndex()[1];
	  WeightZ += it.Get() * (it.GetIndex()[2] + 1);
      sum1    += it.Get();
	 }

	   TotalSum += sum1;
       TotalWeightX += WeightX;
	   TotalWeightY += WeightY;
	   TotalWeightZ += WeightZ;

       p.SetPosition(WeightX/sum1, WeightY/sum1, WeightZ/sum1);
       p.SetRed(1);
       p.SetGreen(0);
       p.SetBlue(0);
       p.SetAlpha(1.0);
       list.push_back(p);
	}
     m_COF[0] = TotalWeightX/TotalSum;
     m_COF[1] = TotalWeightY/TotalSum;
	 m_COF[2] = TotalWeightZ/TotalSum;
     std::cout<<"m_COF:"<<m_COF<<std::endl;
	 m_Blob->GetProperty()->SetName("COFs Blob Object");
     m_Blob->SetId(1);
     m_Blob->SetPoints(list);
    

  filter2->GraftOutput( this->GetOutput() );
  this->GraftOutput(filter2->GetOutput());
    
}


template <class TImageType>
void
FindObjFilter<TImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "Threshold:" << this->m_Threshold
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//FindObjFilter Template Ends
//---------------------------------------------------------------------------------------------------