//#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageSliceConstIteratorWithIndex.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include<math.h> 

#define F1  "Angular Second Moment "
#define F2  "Contrast              "
#define F3  "Correlation           "
#define F4  "Variance              "
#define F5  "Inverse Diff Moment   "
#define F6  "Sum Average           "
#define F7  "Sum Variance          "
#define F8  "Sum Entropy           "
#define F9  "Entropy               "
#define F10 "Difference Variance   "
#define F11 "Difference Entropy    "
#define F12 "Meas of Correlation-1 "
#define F13 "Meas of Correlation-2 "

#define Ng 256

//---------------------------------------------------------------------------------------------------
//CO2Feature3DFilter Template Begins: This filter calculate texture features based on the Co-Occurrence
//matrix
//---------------------------------------------------------------------------------------------------
namespace itk{
template<class TImageType>
class ITK_EXPORT CO2Feature3DFilter : 
  public ImageToImageFilter<TImageType, TImageType>
{
public:
  typedef CO2Feature3DFilter                             Self;
  typedef ImageToImageFilter<TImageType,TImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;
  
  itkNewMacro(Self);
  itkTypeMacro(CO2Feature3DFilter, ImageToImageFilter); 

  void PrintSelf(std::ostream& os, Indent indent) const;
  typedef typename Image<unsigned int, 3> COImageType; 
  typedef typename Image<unsigned int, 2> COImageType1; 
 
  typedef typename COImageType::IndexType COIndexType;
  typedef typename TImageType::PixelType  PixelType;
  typedef typename Matrix<float, Ng, Ng>        MatrixType; // Co-Occurrence Matrix
  typedef typename Vector<float, Ng>       RowColVectorType; //Px(i), Py(i)
  typedef typename Vector<float, 2*Ng + 1>       SumVectorType;    //Px+y(k)   2:2ng (2:512)
  typedef typename Vector<float, Ng>       DiffVectorType;    //Px-y(k)   0:ng-1 (0:255)

  typedef typename Vector<float, 13>      FeatureVectorType;
  typedef typename Vector<float, 26>      TFeatureVectorType;
  typedef typename Matrix<float, 13, 13>      FeatureMatrixType;

  typedef typename NumericTraits< PixelType >::AccumulateType SumType;
  typedef typename NumericTraits< SumType >::RealType MeanType;
  typedef typename TImageType::IndexType IndexType;
  typedef typename TImageType::RegionType RegionType;
  
  typedef typename COImageType::Pointer COType;

  itkGetMacro(CoOccurrence,COType);
  itkSetMacro(CoOccurrence,COType);

  itkGetMacro(RangeVector,FeatureVectorType);
  itkSetMacro(RangeVector,FeatureVectorType);
  
  itkGetMacro(AverageVector,FeatureVectorType);
  itkSetMacro(AverageVector,FeatureVectorType);

  itkGetMacro(FeatureVector,TFeatureVectorType);
  itkSetMacro(FeatureVector,TFeatureVectorType);

protected:

	  CO2Feature3DFilter();

protected:
  typedef ImageRegionIterator<COImageType1> outputIteratorType;
  typedef ImageLinearConstIteratorWithIndex<COImageType1> COIteratorType;
  typedef ImageLinearIteratorWithIndex<COImageType1> COIteratorType1;
  typedef ImageRegionConstIteratorWithIndex<COImageType> IteratorType;
  typedef ImageSliceConstIteratorWithIndex<COImageType> SliceIteratorType;

  void GenerateData();

private:
  CO2Feature3DFilter(Self&);
  void operator=(const Self&);

  COType             m_CoOccurrence;
  FeatureVectorType  m_AverageVector;
  FeatureVectorType  m_RangeVector;
  TFeatureVectorType m_FeatureVector;
  
 };
}

namespace itk
{
 template <class TImageType>
 CO2Feature3DFilter<TImageType>  
 ::CO2Feature3DFilter()
 {

 }

template <class TImageType>
void
CO2Feature3DFilter<TImageType>::
GenerateData()
{  
   FeatureVectorType MaxVector;
   FeatureVectorType MinVector;
//  FeatureVectorType RangeVector;
//   FeatureVectorType AverageVector;
   float max = NumericTraits< float >::max();
   float min = NumericTraits< float >::min();
   MaxVector.Fill(min);
   MinVector.Fill(max);
   m_RangeVector.Fill(0);
   m_AverageVector.Fill(0);

   FeatureMatrixType FeatureMatrix;
   MatrixType       COMatrix; 
   MatrixType       pMatrix;
   RowColVectorType Px;
   RowColVectorType Py;
   SumVectorType    PxPy;
   DiffVectorType   PxMy;
   unsigned long    R = 0;

   PxPy.Fill(0);
   PxMy.Fill(0);


    /////////////////////////////////////////////////////////////
   // Obtain variables that are needed for feature calculation//
  /////////////////////////////////////////////////////////////
   // Px -----------------------------------------------------------------
   
   //  COIteratorType it(m_CoOccurrence, m_CoOccurrence->GetRequestedRegion());	
   SliceIteratorType tdIt(m_CoOccurrence, m_CoOccurrence->GetRequestedRegion()); // silce iterator
   tdIt.SetFirstDirection(0);    
   tdIt.SetSecondDirection(1); // this direction decides the line iteration direction

   int sliceNo = 0; 

for( tdIt.GoToBegin(); !tdIt.IsAtEnd(); tdIt.NextSlice())
 { 
    COImageType1::RegionType region;   // COImageType1: 2D ( 1 of 13 Co-Occurrence matrix)
	COImageType1::IndexType start;
    start[0] = 0;
    start[1] = 0;
    COImageType1::SizeType size;
    size[0]  = 256;
    size[1]  = 256;
    region.SetSize(size);
    region.SetIndex(start);

	// create a 2D image CoOccurrence, which is then used to copy one slice from m_CoOccurrence
	COImageType1::Pointer CoOccurrence = COImageType1::New(); 
    CoOccurrence->SetRegions(region);
    CoOccurrence->Allocate();

    COIteratorType1 outputIt( CoOccurrence ,region);
	outputIt.SetDirection(1);
	outputIt.GoToBegin();
    
	// copy one slice to CoOccurrence
    while( !tdIt.IsAtEndOfSlice() )
	{
	 while( !tdIt.IsAtEndOfLine() )
	 {
	 COImageType::IndexType idx = tdIt.GetIndex();
     COMatrix[idx[0]][idx[1]] = tdIt.Get();
	 outputIt.Set(tdIt.Get());
	 ++outputIt;
	 ++tdIt;
	 }
	 outputIt.NextLine();
	 tdIt.NextLine();
	}

    // iteration begins to calculate texture features based on 1 of 13 Co-Occurrence matrix for 1 of 13 directions
    COIteratorType it(CoOccurrence, CoOccurrence->GetRequestedRegion());
    it.SetDirection(1);

   int i = 0;
   for (it.GoToBegin(); !it.IsAtEnd(); it.NextLine())
   {
	Px[i] = 0;  
    while(!it.IsAtEndOfLine())
	{

	COImageType1::IndexType idx = it.GetIndex();
    COMatrix[idx[0]][idx[1]] = it.Get();

    if(idx[0]==0 || idx[1]==0){++it; continue;}

    R += it.Get();
	Px[i] += it.Get();
	++it;
	}
	i++;
   } 

   Px /= R;
   
       
   // Py -----------------------------------------------------------------
   it.SetDirection(0);
   int ii = 0;
   for (it.GoToBegin(); !it.IsAtEnd(); it.NextLine())
   {
    Py[ii] = 0;
    while(!it.IsAtEndOfLine())
	{
	COImageType1::IndexType idx = it.GetIndex();
	if(idx[0]==0 || idx[1]==0){++it; continue;}

	Py[ii] += it.Get();
	++it;
	}
	ii++;
   }
   Py /= R;

   // ux, uy, sigmax, sigmay, HX, HY ------------------------------------------
   float ux = 0;
   float uy = 0;
   float sigmax = 0;
   float sigmay = 0;
   float HX     = 0;
   float HY     = 0;
   
   for(int i = 1; i<Ng; i++)
   {
   ux += i * Px[i];
   uy += i * Py[i];
   if (Px[i] == 0) { continue; }
   HX += -1 * Px[i] * log(Px[i]);
   }

   for(int i = 1; i<Ng; i++)
   {
  // sigmax += pow(Px[i] - ux,2)/(Ng - 1);
  // sigmay += pow(Py[i] - uy,2)/(Ng - 1);
   sigmax += pow(i - ux,2)/(Ng - 1);
   sigmay += pow(i - uy,2)/(Ng - 1);
   if (Py[i] == 0) { continue; }
   HY += -1 * Py[i] * log(Py[i]);
   }

   sigmax = sqrt(sigmax);
   sigmay = sqrt(sigmay);

   // Pmatrix   -------------------------------------------------------------
   pMatrix  = COMatrix/R;
  
   
   // PxPy, PxMy ------------------------------------------------------------
   for(int i = 1; i< Ng; i++)
   {
	for(int j = 1; j< Ng; j++)
	{
    int k1 = i + j;
	int k2 = abs( i - j);
    PxPy[k1] += pMatrix[i][j]; 
	PxMy[k2] += pMatrix[i][j];
	}
   }
   
   // HXY -----------------------------------------------------------------
   float HXY = 0;

   for(int i = 1; i< Ng; i++)
   {
    for(int j=1; j< Ng; j++)
	{
     if(pMatrix[i][j] == 0)
	 {
	 continue;
	 }
	 HXY += -1 * pMatrix[i][j] * log( pMatrix[i][j] );
	}
   }

   // HXY1 -----------------------------------------------------------------
   float HXY1 = 0;
   for(int i = 1; i< Ng; i++)
   {
    for(int j=1; j< Ng; j++)
	{
     if(pMatrix[i][j] == 0 || Px[i] * Py[j] ==0)
	 {
	 continue;
	 }
	 HXY1 += -1 * pMatrix[i][j] * log( Px[i] * Py[j] );
	}
   }

   // HXY2 -----------------------------------------------------------------
   float HXY2 = 0;
   for(int i = 1; i< Ng; i++)
   {
    for(int j=1; j< Ng; j++)
	{
     if(Px[i] == 0 || Py[j] ==0)
	 {
	 continue;
	 }
	 HXY2 += -1 * Px[i] * Py[j] * log( Px[i] * Py[j] );
	}
   }

    /////////////////////////////////////////////////////////////
   //              13 Textural Features Calculation           //
  /////////////////////////////////////////////////////////////
  // F1  "Angular Second Moment "
  float ASM = 0;
  float correlation = 0;
  float variance = 0;
  float iDiffM   = 0;
  // float entropy  = 0;
   for(int i = 1; i< Ng; i++)
   {
    for(int j = 1; j< Ng; j++)
	{
    ASM         += pow(pMatrix[i][j],2);
	correlation += i * j * pMatrix[i][j];
	variance    += pow( static_cast<double>(i - j), 2) * pMatrix[i][j];
	iDiffM      += pMatrix[i][j] * 1/(1 + pow( static_cast<double>(i - j), 2));
	if (pMatrix[i][j] == 0){ continue; }
  //  entropy     += -1 * pMatrix[i][j] * log( pMatrix[i][j] );
	}
   }

  // F2  "Contrast              "
  float Contrast = 0;
   for(int i = 0 ; i< Ng; i++ )
   {
   Contrast += pow(static_cast<double>(i),2) * PxMy[i];
   }

  // F3  "Correlation           "
  float Correlation = (correlation - ux * uy)/(sigmax * sigmay);

  // F4  "Variance              "
  float Variance = variance;
  
  // F5  "Inverse Diff Moment   "
  float IDiffM   = iDiffM;

  // F6  "Sum Average           "
  // F8  "Sum Entropy           "
  float SumAvg   = 0;
  float SumEnt   = 0;
  for(int i = 2; i < 2 * Ng + 1; i++)
  {
   SumAvg       += i * PxPy[i];  
   if (PxPy[i] == 0){ continue; }
   SumEnt       += -1 * PxPy[i] * log(PxPy[i]);
  }

  // F7  "Sum Variance          "
  float SumVar   = 0;
  for(int i = 2; i < 2 * Ng + 1; i++)
  {
   SumVar       += pow(static_cast<double>(i) - SumEnt, 2) * PxPy[i]; 
  }

  // F9  "Entropy               "
   float Entropy = HXY;

  // F11 "Difference Entropy    "
   float DiffEnt = 0;
   for(int i = 0 ; i< Ng; i++ )
   {
	if (PxMy[i] == 0){ continue; }
   DiffEnt += -1 * PxMy[i] * log(PxMy[i]);
   }

  // F10 "Difference Variance   "
   float DiffVar = 0;
   for(int i = 0 ; i< Ng; i++ )
   {
	if (PxMy[i] == 0){ continue; }
    DiffVar += pow(static_cast<double>(i) - DiffEnt, 2) * PxMy[i];
   }

   // F12 "Meas of Correlation-1 "
   // F13 "Meas of Correlation-2 "
   float maxHXY = HX > HY ? HX:HY;
   float InfoMC1 = (HXY - HXY1) /  maxHXY;
   float InfoMC2 = sqrt(1 - exp(-2 * abs(HXY2 - HXY)));

   FeatureMatrix[sliceNo][0] = ASM;
   FeatureMatrix[sliceNo][1] = Contrast;
   FeatureMatrix[sliceNo][2] = Correlation;
   FeatureMatrix[sliceNo][3] = Variance;
   FeatureMatrix[sliceNo][4] = IDiffM;
   FeatureMatrix[sliceNo][5] = SumAvg;
   FeatureMatrix[sliceNo][6] = SumVar;
   FeatureMatrix[sliceNo][7] = SumEnt;
   FeatureMatrix[sliceNo][8] = Entropy;
   FeatureMatrix[sliceNo][9] = DiffVar;
   FeatureMatrix[sliceNo][10] = DiffEnt;
   FeatureMatrix[sliceNo][11] = InfoMC1;
   FeatureMatrix[sliceNo][12] = InfoMC2;
   

   // update range and average vectors
   for(int i = 0; i<  13; i++)
   {    
	   MaxVector[i]      = FeatureMatrix[sliceNo][i] > MaxVector[i] ? FeatureMatrix[sliceNo][i]:MaxVector[i];
       MinVector[i]      = FeatureMatrix[sliceNo][i] < MinVector[i] ? FeatureMatrix[sliceNo][i]:MinVector[i];
	   m_RangeVector[i]    = MaxVector[i] - MinVector[i];
       m_AverageVector[i] += FeatureMatrix[sliceNo][i] / 13;
   }
   
   sliceNo++;

 }
  
  //std::cout<<"ASM:"<<ASM<<std::endl;
  //std::cout<<"Contrast:"<<Contrast<<std::endl;
  //std::cout<<"Correlation:"<<Correlation<<std::endl;
  //std::cout<<"Variance:"<<Variance<<std::endl;
  //std::cout<<"IDiffM:"<<IDiffM<<std::endl;
  //std::cout<<"SumAvg:"<<SumAvg<<std::endl;
  //std::cout<<"SumVar:"<<SumVar<<std::endl;
  //std::cout<<"SumEnt:"<<SumEnt<<std::endl;
  //std::cout<<"Entropy:"<<Entropy<<std::endl;
  //std::cout<<"DiffVar:"<<DiffVar<<std::endl;
  //std::cout<<"DiffEnt:"<<DiffEnt<<std::endl;
  //std::cout<<"InfoMC1:"<<InfoMC1<<std::endl;
  //std::cout<<"InfoMC2:"<<InfoMC2<<std::endl;
  for(int i = 0; i<  13; i++)
 {
   m_FeatureVector[i] = m_AverageVector[i];
   m_FeatureVector[i + 13] = m_RangeVector[i];
  }

  this->GraftOutput( this->GetOutput() );  
  this->GraftOutput( this->GetOutput() );
}




template <class TImageType>
void
CO2Feature3DFilter<TImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "m_CoOccurrence" << this->m_CoOccurrence
    << std::endl;
}

}

//---------------------------------------------------------------------------------------------------
//CO2Feature3DFilter Template Ends
//---------------------------------------------------------------------------------------------------