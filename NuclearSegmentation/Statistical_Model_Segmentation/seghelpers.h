

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRegionOfInterestImageFilter.h"



void norm_vec_std(double **mat1, long m, long n1, double **mat2, long n2);

int Num_Lines2(char *root, int mode);
FILE *FDeclare2(char *root, char *extension, char key);
void calc_norm(double **mat, long D,long n, double *ave,double *stdev);
std::string convert2string(unsigned long id);
