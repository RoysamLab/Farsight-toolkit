#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include <conio.h> //KREIN SIAD IT BEAKS WITH gcc LINUX
#include <string.h>
#include <time.h>
//#include "itoa.h"
#include <vtkTable.h>
#include <vtkVariant.h>
#include <vtkVariantArray.h>
#include <vtkSmartPointer.h>
#include <vtkVariantArray.h>
#include <vtkDoubleArray.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
//KREIN said to include ita.h for LINUX

#define VERSIONN     0.99
#define LIMIT        5000
#define MEMLIMIT    10000
#define LICENCEE        0
#define LINUX           0
#define NUM_SEN_LEVELS 20
#define SEN_RANGE       2
#define LINUXFLAG       0


#define PW 1302
#define LINUX 0

typedef double *PFLOAT;
typedef PFLOAT VECTOR;
typedef PFLOAT *MATRIX;
typedef MATRIX *TENSOR;

typedef char   *PCHAR;
typedef PCHAR  CVECTOR;
typedef PCHAR  *CMATRIX;

typedef int   *PINT;
typedef PINT  IVECTOR;

typedef long int   *LINT;
typedef LINT  LIVECTOR;

typedef long double *PLFLOAT;
typedef PLFLOAT  DVECTOR;
typedef PLFLOAT  *DMATRIX;





class ClusClus 
{
public:
	ClusClus(vtkSmartPointer<vtkTable> Table);
	~ClusClus();
private :
	int        i, l, k, num_data, num_feats, seed, num_trials, gap_flag;
	int        piv_index, link_mode, num_gaps, naive_plot_flag;
	double     avg, std, off_set, mem, ratio, dispersion;
	VECTOR     disposition;
	MATRIX     my_data, thy_data, metrix, GAP, HASTER, mergers;
	
	vtkSmartPointer<vtkTable> table;
	vtkSmartPointer<vtkTable> ClusterProgressTable;

public:
	FILE	*FDeclare_(char *root, char *extension, char key);
	void	AllocateCols_(PFLOAT matrix[], int nRows, int nCols);
	void	ConvertTxtVtk(std::string hname,vtkSmartPointer<vtkTable> table,int num_row,int num_col);
	void	VectorAllocate_(VECTOR *vector, int nCols);
	void	ConvertTxtVtk(char *root,vtkSmartPointer<vtkTable> table,int num_row,int num_col);
	void	Cluster_Hierarchical4_Master_(char *root); 
	int		Determine_File_Chars__(char *root, int *num_data, int *numfeats);
	int		My_Yes_No_Querry_();
	void	MatrixAllocate_(MATRIX *pmatrix, int nRows, int nCols);
	void	Zero_Matrix_(int N, int M, MATRIX A);
	void	CreateDataMatrix(vtkSmartPointer<vtkTable> table,MATRIX my_data);
	double	Cluster_Hierarchical4_Slave_(int num_data, int num_feats, MATRIX my_data,MATRIX metrix, MATRIX mergers, int link_mode, int print_mode); 
	void	Print_Matrixx_(int num_data, int num_feats, MATRIX gerissa, char *root);
	void	Bolivar_Sort_(int numRows, int numColumns, int sort_column, MATRIX big, int mode);
	void	Generate_Synthetic_Data_(int num_data, int num_feats, int seed,MATRIX my_data, MATRIX thy_data);
	void	Stand_Devv_(int num_data, double *avg, double *std, VECTOR my_vector);
	void	MatrixFree_(MATRIX matrix, int nRows);
	void    Data_Distances_V_(int num_data,int num_feats, MATRIX my_data,VECTOR data_distancer);
	double	Calculate_Proto_Distances_V_(int num_clusters, int num_data, int num_feats,VECTOR num_cluster_data, MATRIX my_data, VECTOR cluster_distancer,VECTOR data_distancer, int link_mode);
	double	Update_Proto_Distances_V_(int num_clusters, int num_data, int num_feats,int pirot1, int pirot2, VECTOR num_cluster_data, MATRIX my_data,VECTOR cluster_distancer, VECTOR data_distancer, int link_mode) ;
	double	Merge_ZE_Clusters_V_(int num_data, int *num_clusters, int num_feats,MATRIX my_data, VECTOR proto_distancer, int *pivot1, int *pivot2);
	double  Distance_(int num_feats, VECTOR A, VECTOR B); 
	
	void	SCommand_(char *string1, char *string2, char *string3, char *string4,char *string5);
	void	Mergers_To_Progress_(int num_data, int num_feats, char *root);
	void	Zero_Vector_(int N, VECTOR A);
	void	GAP1_plot_(int gap_flag);
	void	OptimalLeafOrder1D_(char *root);
	void	GetKids_(int num_data, int k, VECTOR Tnums, VECTOR pickedup, VECTOR kids, MATRIX T); 
	void	Prepare_MATLAB_Trees_(char *root);
	void	Print_IVector_(int num_data, VECTOR y, char *root);
	void	Mergers_2_Chris_(char *root);
	int		FindRowIndex_(int num_data, int col, int c, MATRIX A);
	static double	Minimilian_(double a1, double a2);
	static double	Maximilian_(double a1, double a2);

	
	int		randomm_(int parameter);
	static int		compare_(const void *a, const void *b);
	void	Print_IMatrixx_(int num_data, int num_feats, MATRIX gerissa, char *root);
	void	Read_Meta_Data_(int num_data, int num_feats, MATRIX my_data, char *root,int mode); 
	vtkSmartPointer<vtkTable> GetProgressTable();

};