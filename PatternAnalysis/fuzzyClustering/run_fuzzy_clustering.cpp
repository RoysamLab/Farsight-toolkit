/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#include "fuzzy_clustering.h"
#include <stdlib.h> /* atoi */

int main(int argc, char* argv[])
{
	if(argc<6)
	{
		std::cerr<< "Usage: fuzzy_clustering InputFileName NumberOfClasses ClusterOutputFileName TrainingSetPercentage TrainingOutputFile [options]" << std::endl;	
		return 0;
	}
	
	//Create a fuzzy clustering object
	//for now, use the default options
	FuzzyClustering *FzCl = new FuzzyClustering(atoi(argv[2]));

	//read the input file
	FzCl->ReadProblem(argv[1]);

	//run fuzzy clustering
	FzCl->Run();

	//This is optional and used for validation
	FzCl->Validate();

	//Write clustering output (cluster assignment and membership values)
	FzCl->WriteClusteringOutputToFile(argv[3]);

	//extract training set 
	FzCl->ExtractTrainingSet(atof(argv[4]));

	//write the training set to a file
	FzCl->WriteTrainingSetToFile(argv[5]);

	//No need to delete
	//this is for testing purposes only
	delete FzCl;

	return 1;
}
