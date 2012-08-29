#ifndef ACTIVEVALIDATION_H
#define ACTIVEVALIDATION_H

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vtkSmartPointer.h>
#include <vtkVariant.h>
#include <vtkTable.h>

#include <ClusClus/Kmeans.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <utility>
#include <string>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <boost/random.hpp>

class ActiveValidation
{
private:
	struct strata{
		int Nk;                             //number of element in kth strata
		int nk;                             //number of element sampled in kth strata
		int hk;                             //positive outcomes from nk trials
		int numtobeDrawn;                   //the number to be drawn in current iteration
		int numleft;                        //the number unselected in each strata
		double sigmak;                      //variance of element in kth strata
		double phatk;                       //true mean of drawn elements in kth strata
		double varphatk;                    //variance of phatk
		double fpcvarphatk;                 //the finite population correction factor
		std::vector<int > idsk;             //set of global ids for each strata
		std::vector<int > selectedidsk;     //set of selected local ids for each strata
		std::vector<int > selectingidsk;    //set of local ids to be drawn for each strata
		std::vector<int > unselectedidsk;   //set of unselected local ids for each strata
		vnl_matrix<double> strdata;         //sub data set for each strata
		vnl_vector<int> strlable;           //sub lable set for each strata
	};

public:
	//Constructor
	ActiveValidation();
	//Deconstructor
	~ActiveValidation();
	//Initialize an instance
	void Initializing(vnl_matrix<double> data, vnl_vector<int> lable, int K, double delta, unsigned long seed);
	void Initializing(std::vector<std::vector<double > > datavector, 
		              std::vector<int > lablevector, int K, double delta, unsigned long seed);
	//Get homogenous groups
	void Stratifing();
	//Sampling from original dataset
	void Sampling();
	//Output files
	void FileWrite();

public:
	double phat;                           //true mean of entire dataset
	double varphat;                        //variance of phat
	int numiteration;                      //number of iteration
	int numsampled;                        //total number of samples
	unsigned long seed;                             //seed for random generator

private:
	vnl_matrix<double> data;               //original data set to be processed
	vnl_vector<int> lable;	               //lables to original data
	int N;                                 //total number of elements in the entire dataset                               
	int numfeat;                           //number of features in original data set
	int numbin;                            //number of bins                                                  
	double delta;                          //interval radious for significance level
	std::vector<strata > stratas;          //vector of stratas for entire dataset
                    
	                        

private:
	//Calculate hk for each strata
	void Calculatehk(int kth_strata);
	//Check strata whether it need to to be split
	bool SplitCheck(int kth_strata);
	//Split strata with variance larger than threshold
	void StrataSplit(int kth_strata);
	//Random sampling m variables from n entries
	void RandomSampling(int m, int n, int kth_strata);
	//Calculate sigmak for each strata
	void CalculateSigmak(int kth_strata);
	//Calculate number of elements to be drawn
	void CalculatenumToBeDrawn();
	//Calculate phatk for each strata
	void Calculatephatak(int kth_strata);
	//Calculate phat
	void Calculatephata();
	//Calculate variance of phat
	void CalculateVarphat();
	//Calculate variance of phatk
	void CalculateVarphatk(int kth_strata);
	//Check whether convergence has been derived
	bool ConvergeCheck();
	//Update bins after random sampling from each strata
	void UpdateBins(int kth_strata);
};
#endif