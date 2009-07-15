////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Calculate SLF9 features (DNA )based on Segmentation Result, Translation of Murphy's SLF code "ml_obj2feat.m"///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "itkRescaleIntensityImageFilter.h"

#include "itkBlobSpatialObject.h"

#include "itkNumericTraits.h"
#include <math.h> 

#include "itkVector.h"
#include "itkArray.h"
#include "itkEuclideanDistance.h"
//#include "itkImageRegionConstIterator.h"

//---------------------------------------------------------------------------------------------------
//Obj2FeatureDNAFilter Template Begins
//---------------------------------------------------------------------------------------------------
namespace itk{
template<class TImageType>
class ITK_EXPORT Obj2FeatureDNAFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef Obj2FeatureDNAFilter                               Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro(Obj2FeatureDNAFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;
   
  typedef typename TImageType::PixelType PixelType;

  typedef BlobSpatialObject<3> BlobType;
  typedef BlobType::Pointer BlobPointer;
  typedef SpatialObjectPoint<3> BlobPointType;
  
  typedef typename TImageType::Pointer ImagePointer;

  typedef typename Array< float > MeasurementVectorType;
  typedef typename Vector<float, 28> VectorType;
  typedef typename Vector<float, 3> COFVectorType;
  typedef typename NumericTraits< PixelType >::AccumulateType SumType;

  typedef typename NumericTraits< SumType >::RealType MeanType;
  typedef typename TImageType::IndexType IndexType;
  typedef typename TImageType::RegionType RegionType;
  
  itkGetMacro(BlobDNA,BlobPointer);
  itkSetMacro(BlobDNA,BlobPointer);

  itkGetMacro(TotalVolumeDNA,MeanType);
  itkSetMacro(TotalVolumeDNA,MeanType);

  itkGetMacro(AverageVolumeDNA,MeanType);
  itkSetMacro(AverageVolumeDNA,MeanType);

  itkGetMacro(COFDNA, COFVectorType);
  itkSetMacro(COFDNA, COFVectorType);

  itkGetMacro(NobjectDNA, int);
  itkSetMacro(NobjectDNA, int);
  
  itkGetMacro(FeatureVector,VectorType);
  itkSetMacro(FeatureVector,VectorType);

  itkGetMacro(COFProt,COFVectorType);
  itkSetMacro(COFProt,COFVectorType);
  itkGetMacro(VolProt,float);
  itkSetMacro(VolProt,float);

protected:

	  Obj2FeatureDNAFilter();

protected:

  typedef RescaleIntensityImageFilter< TImageType, TImageType > RescalerType;
  void GenerateData();

private:
  Obj2FeatureDNAFilter(Self&);
  void operator=(const Self&);
  typename RescalerType::Pointer     m_RescaleFilter;

  MeanType               m_AverageVolumeDNA;
  MeanType               m_TotalVolumeDNA;
  typename BlobPointer   m_BlobDNA;    
  COFVectorType          m_COFDNA;
  int                    m_NobjectDNA;
 VectorType              m_FeatureVector;
  COFVectorType          m_COFProt;
  float                  m_VolProt;

 };
}

namespace itk
{
 template <class TImageType>
 Obj2FeatureDNAFilter<TImageType>  
 ::Obj2FeatureDNAFilter()
 {
	 m_RescaleFilter = RescalerType::New();
	 m_RescaleFilter->SetOutputMinimum(0);
     m_RescaleFilter->SetOutputMaximum(255);
 }

template <class TImageType>
void
Obj2FeatureDNAFilter<TImageType>::
GenerateData()
{     
  typedef ImageRegionConstIterator<TImageType> ConstIteratorType;
  float AvgDistance = 0;
  float AvgDistanceXY= 0;
  float AvgDistanceZ= 0;

  // Compute Distance Feature (SLF9_6-8, horizontal-vertical directional SLF9_6-8)
  MeasurementVectorType distanceVector( m_NobjectDNA ); 
  MeasurementVectorType distanceVectorXY( m_NobjectDNA ); 
  MeasurementVectorType distanceVectorZ( m_NobjectDNA ); 
  MeasurementVectorType distanceVectorCOF( 1 ); 
  MeasurementVectorType distanceVectorCOFXY( 1 ); 
  MeasurementVectorType distanceVectorCOFZ( 1 ); 

  MeasurementVectorType queryPoint( 3 );
  MeasurementVectorType queryPointXY( 2 );
  MeasurementVectorType queryPointZ( 1 );
  MeasurementVectorType queryPointCOF( 3 );
  MeasurementVectorType queryPointCOFXY( 2 );
  MeasurementVectorType queryPointCOFZ( 1 );

  typedef Statistics::EuclideanDistance< MeasurementVectorType >
  DistanceMetricType;

  DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
  DistanceMetricType::Pointer distanceMetricXY = DistanceMetricType::New();
  DistanceMetricType::Pointer distanceMetricZ = DistanceMetricType::New();

  DistanceMetricType::OriginType originPoint( 3 );
  DistanceMetricType::OriginType originPointXY( 2 );
  DistanceMetricType::OriginType originPointZ( 1 );

  originPoint[0] = originPointXY[0] = m_COFDNA[0];
  originPoint[1] = originPointXY[1] = m_COFDNA[1];
  originPoint[2] = originPointZ[0] = m_COFDNA[2];

  DistanceMetricType::Pointer distanceMetricCOF = DistanceMetricType::New();
  DistanceMetricType::Pointer distanceMetricCOFXY = DistanceMetricType::New();
  DistanceMetricType::Pointer distanceMetricCOFZ = DistanceMetricType::New();
  distanceMetricCOF->SetOrigin( originPoint );
  distanceMetricCOFXY->SetOrigin(originPointXY);
  distanceMetricCOFZ->SetOrigin(originPointZ);

  queryPointCOF[0] = queryPointCOFXY[0] = m_COFProt[0];
  queryPointCOF[1] = queryPointCOFXY[1] = m_COFProt[1];
  queryPointCOF[2] = queryPointCOFZ[0]  = m_COFProt[2];

  distanceVectorCOF[0]= distanceMetricCOF->Evaluate( queryPointCOF );
  distanceVectorCOFXY[0]= distanceMetricCOFXY->Evaluate( queryPointCOFXY );

  distanceVectorCOFZ[0] = queryPointCOFZ[0] - originPointZ[0];

  BlobType::PointListType::const_iterator Blobit = m_BlobDNA->GetPoints().begin();
  
  int index = 0;
  while(Blobit != m_BlobDNA->GetPoints().end())
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
   AvgDistance = AvgDistance/m_NobjectDNA;
   AvgDistanceXY = AvgDistanceXY/m_NobjectDNA;
   AvgDistanceZ = AvgDistanceZ/m_NobjectDNA;


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

  for(int index1 = 0; index1<= m_NobjectDNA-1; index1++)
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

   SUM += pow( (distanceVector[index1] - AvgDistance), 2)/m_NobjectDNA;
   SUMXY += pow( (distanceVectorXY[index1] - AvgDistanceXY), 2)/m_NobjectDNA;
   SUMZ += pow( (distanceVectorZ[index1] - AvgDistanceZ), 2)/m_NobjectDNA;

  }

  SDDistance = sqrt(static_cast<double>(SUM));
  SDDistanceXY = sqrt(static_cast<double>(SUMXY));
  SDDistanceZ = sqrt(static_cast<double>(SUMZ));

  float Max2MinDistance = maxDistance/minDistance;
  float Max2MinDistanceXY = maxDistanceXY/minDistanceXY;
  float Max2MinDistanceZ = maxDistanceZ/minDistanceZ;
  
  m_FeatureVector[0] =  0;
  m_FeatureVector[1] =  0;
  m_FeatureVector[2] =  0;
  m_FeatureVector[3] =  0;
  m_FeatureVector[4] =  0;
  m_FeatureVector[5] =  0;
  m_FeatureVector[6] =  0;
  m_FeatureVector[7] =  0;
  m_FeatureVector[14] = 0;
  m_FeatureVector[15] = 0;
  m_FeatureVector[16] = 0;
  m_FeatureVector[17] = 0;
  m_FeatureVector[18] = 0;
  m_FeatureVector[19] = 0;

    // DNA features
  m_FeatureVector[8] = AvgDistance;
  m_FeatureVector[9] = SDDistance;
  m_FeatureVector[10] = Max2MinDistanceXY;
  m_FeatureVector[11] = distanceVectorCOF[0]; 
  m_FeatureVector[12] = m_VolProt/m_TotalVolumeDNA;
  m_FeatureVector[13] = 0;            //Fraction of protein fluorescence overlapping with DNA fluorescence, to be calculated at SLF9FeatureCal.cxx
  m_FeatureVector[20] = AvgDistanceXY;
  m_FeatureVector[21] = SDDistanceXY;
  m_FeatureVector[22] = Max2MinDistanceXY;
  m_FeatureVector[23] = AvgDistanceZ;
  m_FeatureVector[24] = SDDistanceZ;
  m_FeatureVector[25] = Max2MinDistanceZ;
  m_FeatureVector[26] = distanceVectorCOFXY[0]; //The horizontal distance between the protein COF and the DNA COF
  m_FeatureVector[27] = distanceVectorCOFZ[0];  //The signed vertical distance between the protein COF and the DNA COF


  m_RescaleFilter->SetInput(this->GetInput());
  m_RescaleFilter->GraftOutput( this->GetOutput() );
  m_RescaleFilter->Update();
  this->GraftOutput(m_RescaleFilter->GetOutput());
    
}


template <class TImageType>
void
Obj2FeatureDNAFilter<TImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "Threshold:" << this->m_NobjectDNA
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//Obj2FeatureDNAFilter Template Ends
//---------------------------------------------------------------------------------------------------
