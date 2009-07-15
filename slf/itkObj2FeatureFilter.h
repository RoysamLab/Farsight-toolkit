//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Calculate SLF9 features based on Segmentation Result, Translation of Murphy's SLF code "ml_obj2feat.m"///
////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "itkBlobSpatialObject.h"
#include "itkNumericTraits.h"
#include <math.h> 

#include "itkVector.h"
#include "itkArray.h"
#include "itkEuclideanDistance.h"
//---------------------------------------------------------------------------------------------------
//Obj2FeatureFilter Template Begins
//---------------------------------------------------------------------------------------------------
namespace itk{
template<class TImageType>
class ITK_EXPORT Obj2FeatureFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef Obj2FeatureFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro(Obj2FeatureFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;
   
  typedef typename TImageType::PixelType PixelType;

  typedef BlobSpatialObject<3> BlobType;
  typedef BlobType::Pointer BlobPointer;
  typedef SpatialObjectPoint<3> BlobPointType;
  
  typedef typename Array< float > MeasurementVectorType;
  typedef typename Vector<float, 28> VectorType;
  typedef typename Vector<float, 3> COFVectorType;
  typedef typename NumericTraits< PixelType >::AccumulateType SumType;

  typedef typename NumericTraits< SumType >::RealType MeanType;
  typedef typename TImageType::IndexType IndexType;

  itkGetMacro(Blob,BlobPointer);
  itkSetMacro(Blob,BlobPointer);

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
  
  itkGetMacro(FeatureVector,VectorType);
  itkSetMacro(FeatureVector,VectorType);
protected:

	  Obj2FeatureFilter();

protected:
  
  void GenerateData();
private:
  Obj2FeatureFilter(Self&);
  void operator=(const Self&);
  int                    m_Threshold;
  MeanType               m_AverageVolume;
  MeanType               m_SDVolume;
  float                  m_Max2Min;
  typename BlobPointer   m_Blob;    
  COFVectorType          m_COF;
  int                    m_Nobject;
 VectorType              m_FeatureVector;
//  TImageType inputImage;

 };
}

namespace itk
{
 template <class TImageType>
 Obj2FeatureFilter<TImageType>  
 ::Obj2FeatureFilter()
 {
 }

template <class TImageType>
void
Obj2FeatureFilter<TImageType>::
GenerateData()
{     
  
  float AvgDistance = 0;
  float AvgDistanceXY= 0;
  float AvgDistanceZ= 0;

  MeasurementVectorType distanceVector( m_Nobject ); 
  MeasurementVectorType distanceVectorXY( m_Nobject ); 
  MeasurementVectorType distanceVectorZ( m_Nobject ); 

  MeasurementVectorType queryPoint( 3 );
  MeasurementVectorType queryPointXY( 2 );
  MeasurementVectorType queryPointZ( 1 );

  typedef Statistics::EuclideanDistance< MeasurementVectorType >
  DistanceMetricType;
  DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
  DistanceMetricType::Pointer distanceMetricXY = DistanceMetricType::New();
  DistanceMetricType::Pointer distanceMetricZ = DistanceMetricType::New();

  DistanceMetricType::OriginType originPoint( 3 );
  DistanceMetricType::OriginType originPointXY( 2 );
  DistanceMetricType::OriginType originPointZ( 1 );

  originPoint[0] = originPointXY[0] = m_COF[0];
  originPoint[1] = originPointXY[1] = m_COF[1];
  originPoint[2] = originPointZ[0] = m_COF[2];
  
  BlobType::PointListType::const_iterator Blobit = m_Blob->GetPoints().begin();

  int index = 0;
  while(Blobit != m_Blob->GetPoints().end())
  {
   queryPoint[0] = queryPointXY[0] = (*Blobit).GetPosition()[0];
   queryPoint[1] = queryPointXY[1] = (*Blobit).GetPosition()[1];
   queryPoint[2] = queryPointZ[0] = (*Blobit).GetPosition()[2];


   distanceMetric->SetOrigin( originPoint );
   distanceMetricXY->SetOrigin( originPointXY );
   distanceMetricZ->SetOrigin( originPointZ );

   distanceVector[index] = distanceMetric->Evaluate( queryPoint );
   distanceVectorXY[index] = distanceMetricXY->Evaluate( queryPointXY );
   distanceVectorZ[index] = distanceMetricZ->Evaluate( queryPointZ );
   
   AvgDistance  += distanceVector[index]; 
   AvgDistanceXY += distanceVectorXY[index]; 
   AvgDistanceZ += distanceVectorZ[index]; 

   index++;
   Blobit++;
  }
   AvgDistance = AvgDistance/m_Nobject;
   AvgDistanceXY = AvgDistanceXY/m_Nobject;
   AvgDistanceZ = AvgDistanceZ/m_Nobject;


  float minDistance = NumericTraits< float >::max();
  float minDistanceXY = NumericTraits< float >::max();
  float minDistanceZ = NumericTraits< float >::max();

  float maxDistance = 0;
  float maxDistanceXY = 0;
  float maxDistanceZ = 0;

  float SDDistance = 0;
  float SDDistanceXY = 0;
  float SDDistanceZ = 0;
  
  float SUM = 0;
  float SUMXY = 0;
  float SUMZ = 0;

  for(int index1 = 0; index1<= m_Nobject-1; index1++)
  {
	  if (distanceVector[index1] != 0){
	  minDistance = minDistance>distanceVector[index1]?distanceVector[index1]:minDistance;
      maxDistance = maxDistance<distanceVector[index1]?distanceVector[index1]:maxDistance;
	  }
	  if (distanceVectorXY[index1] != 0 ){
	  minDistanceXY = minDistanceXY>distanceVectorXY[index1]?distanceVectorXY[index1]:minDistanceXY;
      maxDistanceXY = maxDistanceXY<distanceVectorXY[index1]?distanceVectorXY[index1]:maxDistanceXY;
	  }
	  if (distanceVectorZ[index1] != 0 ){
      minDistanceZ = minDistanceZ>distanceVectorZ[index1]?distanceVectorZ[index1]:minDistanceZ;
      maxDistanceZ = maxDistanceZ<distanceVectorZ[index1]?distanceVectorZ[index1]:maxDistanceZ;
	  }

   SUM += pow( (distanceVector[index1] - AvgDistance), 2)/m_Nobject;
   SUMXY += pow( (distanceVectorXY[index1] - AvgDistanceXY), 2)/m_Nobject;
   SUMZ += pow( (distanceVectorZ[index1] - AvgDistanceZ), 2)/m_Nobject;

  }

  SDDistance = sqrt(static_cast<double>(SUM));
  SDDistanceXY = sqrt(static_cast<double>(SUMXY));
  SDDistanceZ = sqrt(static_cast<double>(SUMZ));

  float Max2MinDistance = maxDistance/minDistance;
  float Max2MinDistanceXY = maxDistanceXY/minDistanceXY;
  float Max2MinDistanceZ = maxDistanceZ/minDistanceZ;


  m_FeatureVector[0] =  m_Nobject;
  m_FeatureVector[1] =  m_Nobject;          // Find a way to calculate holes in the Image
  m_FeatureVector[2] =  m_AverageVolume;
  m_FeatureVector[3] =  m_SDVolume;
  m_FeatureVector[4] =  m_Max2Min;
  m_FeatureVector[5] =  AvgDistance;
  m_FeatureVector[6] =  SDDistance;
  m_FeatureVector[7] =  Max2MinDistance;
  m_FeatureVector[14] = AvgDistanceXY;
  m_FeatureVector[15] = SDDistanceXY;
  m_FeatureVector[16] = Max2MinDistanceXY;
  m_FeatureVector[17] = AvgDistanceZ;
  m_FeatureVector[18] = SDDistanceZ;
  m_FeatureVector[19] = Max2MinDistanceZ;

    // DNA features
  m_FeatureVector[8] = 0;
  m_FeatureVector[9] = 0;
  m_FeatureVector[10] = 0;
  m_FeatureVector[11] = 0;
  m_FeatureVector[12] = 0;
  m_FeatureVector[13] = 0;
  m_FeatureVector[20] = 0;
  m_FeatureVector[21] = 0;
  m_FeatureVector[22] = 0;
  m_FeatureVector[23] = 0;
  m_FeatureVector[24] = 0;
  m_FeatureVector[25] = 0;
  m_FeatureVector[26] = 0;
  m_FeatureVector[27] = 0;

  this->GraftOutput(this->GetOutput());
    
}


template <class TImageType>
void
Obj2FeatureFilter<TImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "Threshold:" << this->m_Threshold
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//Obj2FeatureFilter Template Ends
//---------------------------------------------------------------------------------------------------
