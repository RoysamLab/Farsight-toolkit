#include "fuzzy_clustering.h"

int main(int argc, char* argv[])
{
	if(argc<6)
	{
		cout<< "Usage: fuzzy_clustering InputFileName NumberOfClasses ClusterOutputFileName TrainingSetPercentage TrainingOutputFile [options]" << 
endl;	
		return 0;
	}
	
	//Create a fuzzy clustering object
	//for now, use the default options
	FuzzyClustering *FzCl = new FuzzyClustering(atoi(argv[2]));

	//read the input file
	FzCl->ReadProblem(argv[1]);

	//run fuzzy clustering
	FzCl->Run();

	//Write clustering output (cluster assignment and membership values)
	FzCl->WriteClusteringOutputToFile(argv[3]);

	//extract training set 
	FzCl->ExtractTrainingSet(atoi(argv[4]));

	//write the training set to a file
	FzCl->WriteTrainingSetToFile(argv[5]);

	//No need to delete
	//this is for testing purposes only
	delete FzCl;

	return 1;
}