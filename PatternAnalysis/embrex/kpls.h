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

#ifndef KPLS_H
#define KPLS_H

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifdef USE_KPLS
#include <vcl_vector.h>

typedef double *EMB_PFLOAT;
typedef EMB_PFLOAT VECTOR;
typedef EMB_PFLOAT *MATRIX;

class KPLS
{
public:
	KPLS();
	~KPLS();

	//USE THES FUNCTIONS TO MANUALLY PREPARE THE DATA:
	MATRIX GetDataPtr(int rows, int columns);
	VECTOR GetIDPtr(void);
	VECTOR GetTrainingPtr(void);
	void InitVariables(void);

	MATRIX myData_dup;
	int num_classes;	//Number of Classes
	VECTOR myKnownClass;	//If we load data for which we know the class it goes here

	//USE THESE FUNCTIONS TO LOAD FROM A FILE
	void LoadFromMetaNeuralFormat(char * train_fname, char * test_fname);
    void LoadData(char * train_fname,int sz);
	void LoadNodeFeatures(std::vector<double> in_features );
	void LoadAllFeatures(std::vector< std::vector< double > > infeatures );


	//SETUP PARAMETERS:
	void SetLatentVars(int n) { latent_vars = n; };
	void SetSigma(int n) { sigma = n; };

	//PROCESS:	
	void ScaleData();
	void ScaleData2(double ** myData);
	void Train();
	void Classify();
	void Statistics();
	void SaveResults();
	void SaveLatents();
	std::vector< std::vector<double> > Classify2();

	//GET RESULTS:
	VECTOR GetPredictions(void);

	struct {bool scaled; bool trained; bool classified;} flag;

private:
	void Free_All(void);
	void Free_Intermediates(void);
	void Allocate_Intermediates(void);
	void ScaleMatrix(MATRIX data, MATRIX scalers, int rows, int columns);
	void DeScaleMatrix(MATRIX data, MATRIX scalers, int rows, int columns);
	//void Predict(VECTOR ids, MATRIX targets, MATRIX predictions, int size);
	void Predict(MATRIX predictions);
	std::vector< std::vector<double> > Predict2(MATRIX predictions);

	int latent_vars;	//default to 5
	int sigma;			//default to 20

	//PRIVATE
	VECTOR myID;			//We associate with an ID
	MATRIX myData;			//All Data ( scaled in ScaleData, if called )
	VECTOR myTraining;		//Training Class membership (-1 for not-assigned)
	MATRIX myResponses;		//Provides Responses for each object to each class.
	VECTOR myPredictions;	//Predicted class membership
	//VECTOR myKnownClass;	//If we load data for which we know the class it goes here
							//This is used for Checking our classification for accuracy
	MATRIX myFeatScalers;	//HOLDS avg and std devs of each feature

	int num_data;		//Number of Data points
	int num_feats;      //Number of Features (does not count IDs or Class)
	//int num_classes;	//Number of Classes
	int num_train_data; //Number of Data used in training

	//INTERMEDIATE VARIABLES (variables I need to keep around for classification):
	VECTOR tID;					//IDs of data with known classes;
	MATRIX tData;				//DATA WITH KNOWN CLASSES (RESPONSES) FOR TRAINING
	VECTOR ccmatrixx;			//kernel: centering factor for each row
	MATRIX bbmatrixx;			//kernel: weights for each class for each row
	MATRIX tClassScalers;		//HOLDS avg and std devs of each class
};

#endif

#endif
