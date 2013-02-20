/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#include "WeightedMahalsnobisDistance.h"

//int main(int argc, char **argv)
int MDLClassifier::MeanVectorandVarianceMatrix(char *filename)
{
 
  FILE *SampleSpineFile;
  std::string infilename1;
  infilename1 =  filename;
  if(( SampleSpineFile=fopen(infilename1.c_str(), "rb")) == NULL)  // open RealSpine file
    {
    cerr << "couldn't open volume file " << infilename1 << " for input" << endl;
    return -1;
    }
   char str[200];
   int i,j;
   for (i=0;i<FeatureNumber+1;i++) // read the first line;
   {
   if( fscanf(SampleSpineFile,"%s",str) == EOF )
    {
    cerr << "fscanf encountered end of file!" << endl;
    }
   }
   
  typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
  SampleType::Pointer sample = SampleType::New();
  sample->SetMeasurementVectorSize( FeatureNumber );
 

  float temp;
  int IDNumber =0;
  int TotalNumber =0;
  
  while (1)
  {
  TotalNumber ++;
 
  if( fscanf (SampleSpineFile,"%d",&IDNumber) == EOF )
     {
     cout << "The Sample Spine file for Machine Learning is read over!!" << endl;
	 break;
     }
  
  MeasurementVectorType mv;
  
  for (i=0;i<FeatureNumber;i++)
  {
	 if( fscanf (SampleSpineFile,"%f",&temp) == EOF )
     {
     cout << "The Sample Spine file for Machine Learning is read over!" << endl;
	 break;
     }
	 else
	{
	 mv[i] = temp;
	 sample->PushBack( mv ); 
	}
  }
  } // end while 

  fclose (SampleSpineFile);

  typedef itk::Statistics::WeightedMeanCalculator< SampleType >
    WeightedMeanAlgorithmType;
  


  WeightedMeanAlgorithmType::WeightArrayType weightArray( sample->Size() );
 
  weightArray.Fill( 0.5 );
  weightArray[2] = 0.01;
  weightArray[4] = 0.01;

  WeightedMeanAlgorithmType::Pointer weightedMeanAlgorithm = 
                                              WeightedMeanAlgorithmType::New();

  weightedMeanAlgorithm->SetInputSample( sample );
  weightedMeanAlgorithm->SetWeights( &weightArray );
  weightedMeanAlgorithm->Update();

  
  typedef itk::Statistics::WeightedCovarianceCalculator< SampleType >
                                              WeightedCovarianceAlgorithmType;
  
  WeightedCovarianceAlgorithmType::Pointer weightedCovarianceAlgorithm = 
                                        WeightedCovarianceAlgorithmType::New();

  weightedCovarianceAlgorithm->SetInputSample( sample );
  weightedCovarianceAlgorithm->SetMean( weightedMeanAlgorithm->GetOutput() );
  weightedCovarianceAlgorithm->SetWeights( &weightArray );
  weightedCovarianceAlgorithm->Update();

  
  ExampleWeightFunction::Pointer weightFunction = ExampleWeightFunction::New();

  weightedMeanAlgorithm->SetWeightFunction( weightFunction );
  weightedMeanAlgorithm->Update();

  weightedCovarianceAlgorithm->SetMean( weightedMeanAlgorithm->GetOutput() );
  weightedCovarianceAlgorithm->SetWeightFunction( weightFunction );
  weightedCovarianceAlgorithm->Update();
   

  weightedCovarianceAlgorithm->SetMean( 0 );
  weightedCovarianceAlgorithm->SetWeightFunction( weightFunction );
  weightedCovarianceAlgorithm->Update();
  itk::Statistics::WeightedCovarianceCalculator<SampleType>::MeanType* weighted_mean;
  weighted_mean = weightedCovarianceAlgorithm->GetMean();
  
  std::cout << " The learning Feature Weighted-Mean Discriminal Vector is: " << std::endl;
  for(i=0;i<FeatureNumber; i++)
  {
	  MeanVector[i] = weighted_mean->GetElement(i);
	  
  }

  std::cout << std::endl;
  std::cout << " The learning Feature Weighted-Covariance Discriminal Vector is: " ;

  itk::Statistics::WeightedCovarianceCalculator<SampleType>::OutputType* weighted_Covariance;
  weighted_Covariance = weightedCovarianceAlgorithm->GetOutput();
 
 
 

  int IDX =0;
  for(i=0;i<FeatureNumber;i++)
  {
	 for (j=0;j<FeatureNumber; j++)
	 {
     IDX= i*FeatureNumber + j;
	 //InverseCovarianceMatrix[i][j] = weighted_Covariance->GetInverse().get(i,j);
	 *(InverseCovarianceMatrix+IDX) = weighted_Covariance->GetInverse().get(i,j);
	 //std::cout << *(InverseCovarianceMatrix+IDX) << "  ";
	 }
   }
  return 1;
}

double  MDLClassifier::MahalanobisDist(double *sample)
{ 
   double mahalanobis_dist =0;
   int i,j;
   double *DifferenceVector;
   double *VevectorMultiplyCovarianceMatrix;
   DifferenceVector = new double [this->FeatureNumber];
   VevectorMultiplyCovarianceMatrix = new double [this->FeatureNumber];
   
   for (i=0;i<this->FeatureNumber;i++)
   {DifferenceVector[i] = (sample[i]-MeanVector[i]);
    VevectorMultiplyCovarianceMatrix[i] = 0.; 
	// std::cout << DifferenceVector[i] <<std::endl;
   }
   
   int IDX; 
   for (i=0;i<this->FeatureNumber;i++)
   {
     for (j=0;j<this->FeatureNumber;j++)   
     { IDX =  i * this->FeatureNumber  + j;
	   VevectorMultiplyCovarianceMatrix[i] += DifferenceVector[i]*this->InverseCovarianceMatrix[IDX];
	   //std::cout << *(InverseCovarianceMatrix+IDX) << "  ";
     }

	 //std::cout << std::endl << VevectorMultiplyCovarianceMatrix[i] << "  " ;
   };
    
   for (i=0;i<this->FeatureNumber;i++)
   {
     mahalanobis_dist += VevectorMultiplyCovarianceMatrix[i]*DifferenceVector[i];
   }
   
   delete []DifferenceVector;
   delete []VevectorMultiplyCovarianceMatrix;
   return mahalanobis_dist;
}


double  MDLClassifier::MahalanobisDist(double meanDensityBranch, double length_leaf, double meanVesselBranch, int spineOne) 
{
  double x1, x2, x3;
  double mahalanobis_dist = 0;
  /* For dataset time330
  if (spineOne == 1)  {
        x1 = meanDensityBranch - 58.4550;  // minus mean_feature from matlab 
        x2 = length_leaf -    14.0834;              //if (x2>10)  x2=10; //keep long dendrite 
        x3 = meanVesselBranch -  205.4605;
    mahalanobis_dist = x1*x1*0.0032+ 2*x1*x2*(0.0018)+ 2*x1*x3*(-0.001)+ x2*x2*0.0562+ 2*x2*x3*(-0.0017) +x3*x3*0.0006;
  }
  else  {
        x1 = meanDensityBranch - 13.1488;  // assume non-spine has close-to-zero distribution, but the same variance as spines 
        x2 = length_leaf -     6.9848;                             
        x3 = meanVesselBranch -   59.0966;
    mahalanobis_dist = x1*x1*0.0101+ 2*x1*x2*(0.0011)+ 2*x1*x3*(-0.0025)+ x2*x2*0.0202+ 2*x2*x3*(-0.0011) +x3*x3*0.0009;
  }
  */

  // For dataset Trach6A

 // For dataset Trach6A
	/*
	if (spineOne == 1)  {
        x1 = meanDensityBranch - 44.54;  // minus mean_feature from matlab 
        x2 = length_leaf - 14.43;             	//if (x2>10)  x2=10; //keep long dendrite 
        x3 = meanVesselBranch - 198.62;
		mahalanobis_dist = x1*x1*0.0019+ 2*x1*x2*(0.0031)+ 2*x1*x3*(-0.0002)+ x2*x2*0.0252+ 2*x2*x3*0.0002	+x3*x3*0.0004;
	}
	else	{
        x1 = meanDensityBranch - 37.52;  // assume non-spine has close-to-zero distribution, but the same variance as spines 
        x2 = length_leaf - 9.42;                             
        x3 = meanVesselBranch - 178.26;
		mahalanobis_dist = x1*x1*0.0018+ 2*x1*x2*(-0.0002)+ 2*x1*x3*(-0.0003)+ x2*x2*0.0182+ 2*x2*x3*0.0010	+x3*x3*0.0004;
	}
	
*/
  
  // For dataset MBFsp

  if (spineOne == 1)  {
        x1 = meanDensityBranch - 58.4550;  // minus mean_feature from matlab 
        x2 = length_leaf -    14.0834;              //if (x2>10)  x2=10; //keep long dendrite 
        x3 = meanVesselBranch -  58.4550;
    mahalanobis_dist = x1*x1*0.0032+ 2*x1*x2*(0.0018)+ 2*x1*x3*(-0.001)+ x2*x2*0.0562+ 2*x2*x3*(-0.0017) +x3*x3*0.0006;
  }
  else  {
        x1 = meanDensityBranch - 13.1488;  // assume non-spine has close-to-zero distribution, but the same variance as spines 
        x2 = length_leaf -     6.9848;                             
        x3 = meanVesselBranch -   13.1488;
    mahalanobis_dist = x1*x1*0.0101+ 2*x1*x2*(0.0011)+ 2*x1*x3*(-0.0025)+ x2*x2*0.0202+ 2*x2*x3*(-0.0011) +x3*x3*0.0009;
  }
  
  return mahalanobis_dist;
}

/*
int main(int argc, char **argv)
{

 if (argc < 2)
    {
       cerr << "Usage: " << argv[0] << " <RealSpineFile> <NonSpineFile> " << endl;
	   return 1; 
    }

 MDLClassifier Test(3);

 int t= Test.MeanVectorandVarianceMatrix(argv[1]);
 Test.print();

 double *sample;
 sample = new double [3];
 sample[0] = 124; sample[1] = 4.32;  sample[2] = 145.34; 
 double tmp;
 tmp = Test.MahalanobisDist(sample);
 std::cout << tmp;

 getchar();
}
*/